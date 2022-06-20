//
//  main.cpp
//  ColGenFS
//
//  Created by StephanP on 10/21/21.
//  Copyright © 2021 Stephan Patterson. All rights reserved.
//
//  CG on LP(general) copied from ColGenDelFS
//  Version for using the new pricing LP
//

#include "/Library/gurobi801/mac64/include/gurobi_c++.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include "SuppPt.hpp"
#include <numeric>

bool validateDist(const std::vector<SuppPt>::iterator, const std::vector<SuppPt>::iterator);
void mvmult(const int , const int long , const int long , const std::vector< double> &, double * const, const int * const , const int );
void mvmult(const int , const int long , const int long , const std::vector< double> &, std::vector<int> &, double * const, const int * const , const int );
std::vector< double> mvmult(std::vector< std::vector<int> >::iterator,const int, const std::vector< double> &);
std::vector< double> vmmult(const std::vector< std::vector<int> > & , const int , const long int , const double* const );
int vmmult(const std::vector< std::vector<int> > &, const int , const long int , const double* const , std::vector<double> &);
std::vector<double> dualtotals(const double *const, const int, const long int, const int *, const int );
int adjusttotals(std::vector<double> &, const double * const, const int, const long int, const int * const);
int makegreedy(std::vector<SuppPt> &, std::vector<int> & , int *, const int &, std::vector<double> &, const long int &);

int main(int argc, const char * argv[]) {
   
   auto start = std::chrono::steady_clock::now();
   auto t = start;
   auto t0 = t;
   std::vector< std::time_t > pricingtimes;

   std::cout.precision(10);
   int Pnum =4;//If NumMoMeas = 1, Number of months to process. Will go consecutively from start month, skipping any empty months.
   std::string startmonth = "Jul";
   int startyear = 15;

   int NumMoMeas = 1;//Number of Months in a single Measure. Use 1 for original approach of each month is its own measure
   //Otherwise, Number of measures that will be created. Will combine NumMoMeas of months into each measure
//   int fracsol = 0;
   srand(time(NULL));
   
   int warmstart = 1; //If 1, the previous solution is saved and reintroduced before ->optimize()
   // This should be 1 for fastest running times
   int pwarmstart = 0;
   
   int Psnum[Pnum];
   std::vector<SuppPt> Psupp;
   //For the new LP, it makes more sense to have a two-indexed vector PsuppT instead of one long vector Psupp. Converting fully to PsuppT has not yet been implemented 12/15/21
   std::vector< std::vector<SuppPt> > PsuppT(Pnum, std::vector< SuppPt > ((0)));
   int totalsupp = 0;
   int startindices[Pnum];
   int indices[Pnum];
   int endindices[Pnum];
   double lambda[Pnum];

   std::string temp;

   std::ifstream indata;
   std::ostringstream fname;
   int year;
   std::string month;
   fname << "/Users/spatterson/ForXcode/Barycenter/DenverCrime2.csv";//File 2 has months with at least 2 murders only
   indata.open(fname.str());
   if (!indata.good())
   {
      std::cout << "File not found." << std::endl;
      return 1;
   }
   getline(indata, temp, '\n');// Advance past headers
   indata >> year;
   if (startyear < year)
   {
      std::cout << "Invalid Year Selection." << std::endl;
      return 1;
   }
   while (year < startyear)
   {
      getline(indata, temp, '\n');// Advance to next entry
      indata >> year;
   }
   indata >> month;
   while (month != startmonth)
   {
      getline(indata, temp, '\n');
      indata >> year;
      indata >> month;
      if (year != startyear)
      {
         std::cout << "Selected Month/Year not in data." << std::endl;
         return 1;
      }
      if (indata.eof())
      {
         indata.close();
         std::cout << "Invalid Input. End of file reached." << std::endl;
         return 1;
      }
   }
   
   double loc1;
   double loc2;
   indata >> loc1;
   indata >> loc2;
   std::string currentmonth = startmonth;
         
   int tempind = 0;
   for (int i = 0; i < Pnum; ++i)
   {
      Psnum[i] = 0;
      startindices[i] = totalsupp;
      indices[i] = totalsupp;
      for (int j = 0; j < NumMoMeas; ++j)
      {
         while (month == currentmonth )
         {
            //adding shift so that all coordinates are positive
            Psupp.push_back(SuppPt(loc1+115, loc2 , 1.0, totalsupp));
            PsuppT[i].push_back(SuppPt(loc1+115, loc2, 1.0,totalsupp));
            ++tempind;
            ++Psnum[i];
            ++totalsupp;
                  
            indata >> year;
            indata >> month;
            indata >> loc1;
            indata >> loc2;
            if (indata.eof())
            {
               indata.close();
               std::cout << "Invalid Input. End of file reached." << std::endl;
               return 1;
            }
         }
         currentmonth = month;
      }
      std::cout << "Month measure: " << i+1 << " Size: " << Psnum[i] << std::endl;
      double totalmass = 0;
      for (int j = 0; j < Psnum[i]-1; ++j)
      {
         Psupp[startindices[i]+j].mass /= Psnum[i];
         totalmass += Psupp[startindices[i]+j].mass;
      }
      Psupp[totalsupp-1].mass = 1-totalmass;
      currentmonth = month;
      endindices[i] = totalsupp-1;
   }
   Psupp.resize(totalsupp);
   indata.close();
   
   double sum = 0;
   for (int i = 0; i < Pnum-1; ++i)
   {
      lambda[i] = (double)Psnum[i]/totalsupp;
      sum += lambda[i];
   }
   lambda[Pnum-1] = 1-sum;

   
   int index = 0;
   // ensure P_i sum to exactly 1
   std::vector<SuppPt>::iterator Psupptest = Psupp.begin();
   for (int i = 0; i < Pnum; ++i)
   {
      validateDist(Psupptest, Psupptest+Psnum[i]);
      Psupptest += Psnum[i];
   }
   
   std::vector<SuppPt>::iterator Psuppit = Psupp.begin();
/*   for (int j = 0; j < totalsupp; ++j)
   {
      Psuppit->combos = j;
      ++Psuppit;
   }
   Psuppit = Psupp.begin();
   
   index= 0;
   for (int i = 0; i < Pnum; ++i)
   {
      std::cout << "P: " <<i <<std::endl;
      for (int k = 0; k < Psnum[i]; ++k)
      {
         std::cout << Psupp[index].combos <<std::endl;
         ++index;
      }
   }*/

   //Calculate number of possible supp pts S0
   long int S0 = 1;
   for (int i = 0; i < Pnum; ++i)
   {
      S0 *= Psnum[i];
   }
   std::cout << "Size of S0: " << S0 << std::endl;
   
   //Copying Mass and Index information into the (new) 2-indexed support vector
   //Check if this is still needed 5/18/22
   tempind = 0;
   for (int i = 0; i < Pnum; ++i)
   {
      for (int j = 0; j < Psnum[i]; ++j)
      {
         PsuppT[i][j].mass = Psupp[tempind].mass;
         PsuppT[i][j].combos = Psupp[tempind].combos;
//         std::cout << PsuppT[i][j] << " " << Psupp[tempind] <<std::endl;
         ++tempind;
      }
   }
   
   //vector for storing costs of combinations. Very large.
   std::vector<double> c(S0,0);

   //Model for Master Problem
   GRBEnv* env = new GRBEnv();
   GRBModel* model = new GRBModel(*env);
   model->set(GRB_IntParam_Method, 0);
   model->set(GRB_IntParam_Presolve, 0);
   model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
//   model->set(GRB_IntParam_OutputFlag, 0); //Turning off display on screen. Disable to see iterations, objective value, etc.
   
   std::vector< GRBVar > w;
   GRBLinExpr  exp[totalsupp];
   
   index = 0;
   
   for (unsigned long int j = 0; j < S0; ++j)
   {
      double sum1 = 0;
      double sum2 = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         index = indices[i];
         sum1 += lambda[i]*Psupp[index].loc1;
         sum2 += lambda[i]*Psupp[index].loc2;
      }
         
      SuppPt Pbar0 = SuppPt(sum1, sum2, 0.0);
      for (int i = 0; i < Pnum; ++i)
      {
         index = indices[i];
         c[j] += lambda[i]*((Pbar0.loc1-Psupp[index].loc1)*(Pbar0.loc1-Psupp[index].loc1) +(Pbar0.loc2-Psupp[index].loc2)*(Pbar0.loc2-Psupp[index].loc2));
      }
         
      int k = Pnum-1;
      if (indices[k] < endindices[k])
      {
         ++indices[k];
      }
      else
      {
         int temp = k-1;
         while (temp >= 0)
         {
            if (indices[temp] == endindices[temp])
            {
               temp -= 1;
            }
            else
            {
               ++indices[temp];
               break;
            }
         }
         for (int l = k; l > temp; --l)
         {
            indices[l] = startindices[l];
         }
      }
   }
   
   for (int i = 0; i < Pnum; ++i)
   {
      indices[i] = startindices[i];
   }
   
   //For ensuring a column does not get added repeatedly
   std::vector<bool> colin(S0,false);
   //Temporary copy of the support vector, will modify it as mass is removed during greedy construction
   std::vector<SuppPt> Psupp2(totalsupp);
   std::copy(Psuppit, Psuppit+totalsupp, Psupp2.begin() );
   //this method of winit is very memory inefficient
   std::vector<double> winit(S0,0);
   while (Psupp2[Psnum[0]-1].mass >1e-15)
   {
      double minmass = 1;
      for (int i = 0; i < Pnum; ++i)
      {
         if (Psupp2[indices[i]].mass < minmass)
         {
            minmass = Psupp2[indices[i]].mass;
         }
      }
      std::cout << "Minimum mass is " << minmass <<std::endl;
      
      long int loocur = S0;
      int unki = 0;
      int index = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         unsigned long int index2 = Psupp2[indices[i]].combos-index;
         std::cout << Psupp2[indices[i]].combos << " " << index << std::endl;
         std::cout << "Index in Pi: " << Psupp2[indices[i]].combos - index <<std::endl;
         loocur /= Psnum[i];
         unki += loocur*index2;
         index += Psnum[i];
         Psupp2[indices[i]].mass -= minmass;
         if (Psupp2[indices[i]].mass <= 1e-15)
         {
            ++indices[i];
         }
      }
      for (int i = 0; i < Pnum; ++i)
      {
         std:: cout << Psupp2[indices[i]].mass << " ";
      }
      std::cout << std::endl;
      
      std::cout << "current index is " << unki << " out of " << S0-1 <<std::endl;
      winit[unki] = minmass;
   }
   
   //Could now free up Psupp2
   //Begin: build master problem from greedy solution
   long int k = 0;
   int jindex = 0;
   long int Pprod;
   std::vector<int> vbases;
   double oldobj;
   double newobj;
    
   for (std::vector<double>::iterator wit = winit.begin(); wit != winit.end(); ++wit)
   {
         if (*wit != 0)
         {
            w.push_back(model->addVar(0.0, GRB_INFINITY, c[k], GRB_CONTINUOUS));//This updates the objective
            vbases.push_back(0);
            colin[k] = true;
            Pprod = S0/Psnum[0];
            jindex = floor(k/Pprod);
            exp[jindex] += *(w.end()-1);
            for (int l = 1; l < Pnum-1; ++l)
            {
               Pprod /= Psnum[l];
               jindex = startindices[l]+floor( (k % (Pprod*Psnum[l]))/Pprod);
               exp[jindex] += *(w.end()-1);
            }
            jindex = startindices[Pnum-1]+ k % Psnum[Pnum-1];
            exp[jindex] += *(w.end()-1);
         }
      ++k;
   }
   for (int j = 0; j < totalsupp; ++j)
   {
      model->addConstr(exp[j] == Psupp[j].mass);
   }
   winit.resize(0);
   
   t = std::chrono::steady_clock::now();
   std::cout << "Setup Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-start).count()  <<"ms" << std::endl;
   std::cout << "Initial Number of Columns: " << w.size() <<std::endl;
   int iter = 1;
   
   model->optimize(); //First solve for dual values
   newobj = model->get(GRB_DoubleAttr_ObjVal);
   oldobj = newobj;
   t0 = t;
   t = std::chrono::steady_clock::now();
   std::cout << "First LP Solve Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-t0).count()  <<"ms" << std::endl;
 
   int cbases[totalsupp];
   GRBConstr* constrs =model->getConstrs();
   if (warmstart == 1)
   {
      for (int i = 0; i < w.size(); ++i)
      {
         vbases[i] = (w[i].get(GRB_IntAttr_VBasis));
      }
      for (int i = 0; i < totalsupp; ++i)
      {
         cbases[i] = constrs[i].get(GRB_IntAttr_CBasis);
      }
   }
   
   std::vector<int> vbstart = vbases;
   //Beginning new pricing LP here
   //Create a new model for new pricing problem
   GRBModel* modelprice = new GRBModel(*env);
   modelprice->set(GRB_IntParam_Method, 0);
   modelprice->set(GRB_IntAttr_ModelSense, -1);//set to maximize
   modelprice->set(GRB_IntParam_OutputFlag, 0); //Turning off display on screen. Disable to see iterations, objective value, etc.
   
   //Solver Options
   //Set these three
   modelprice->set(GRB_DoubleParam_Heuristics,0); //Turns off heuristics
   modelprice->set(GRB_IntParam_Cuts,0);
   modelprice->set(GRB_IntParam_MIPFocus,2); //Tells the solver to focus on proving optimality over finding feasible solutions
   //These three make little difference, leaving enabled
//   modelprice->set(GRB_IntParam_RINS, 0);
//   modelprice->set(GRB_IntParam_Symmetry, 0);
//   modelprice->set(GRB_IntParam_Disconnected, 0);
   //Leave these on, definitely worse without
//   modelprice->set(GRB_IntParam_Aggregate,0);
   //   modelprice->set(GRB_IntParam_Presolve, 0);

   //Create two sets of variables: 2-index and 4-index
   int Psmax = *std::max_element(Psnum, Psnum+Pnum);
   std::vector< std::vector< GRBVar >> z2(Pnum, std::vector< GRBVar > ((Psmax)));
   std::vector< std::vector< std::vector< std::vector< GRBVar >> >>z4(Pnum,std::vector< std::vector< std::vector< GRBVar >>>(Psmax, std::vector< std::vector< GRBVar >>(Pnum, std::vector< GRBVar >(Psmax))));
   std::vector<int> vbasesp;

   //Create new constraints
   GRBLinExpr expp2[totalsupp*totalsupp*2];
   
   //These vectors are the coefficients for the middle term of the objective
   //The last terms of zik and the first terms of zjm will be zero
   std::vector< std::vector< double > > zikcost(Pnum, std::vector<double> ((Psmax),0));
   std::vector< std::vector< double > > zjmcost(Pnum, std::vector<double> ((Psmax),0));
   for (int i = 0; i < Pnum-1; ++i)
   {
      for (int j = i+1; j < Pnum; ++j)
      {
         for (int k = 0; k < Psnum[i]; ++k)
         {
            double prod = PsuppT[i][k]*PsuppT[i][k];
            zikcost[i][k] += lambda[i]*lambda[j]*prod;
         }
      }
   }
   for (int i = 0; i< Pnum-1; ++i)
   {
      for (int j = i+1; j < Pnum; ++j)
      {
            for (int m = 0; m < Psnum[j]; ++m)
            {
               double prod = PsuppT[j][m]*PsuppT[j][m];
               zjmcost[j][m] += lambda[i]*lambda[j]*prod;
            }
      }
   }
   //Add two-indexed variables to the LP
   int currentconst = 0;
   int totalconstraints = 0;
   for (int i = 0; i < Pnum; ++i)
   {
      GRBLinExpr expp;
      for (int j = 0; j < Psnum[i]; ++j)
      {
         std::stringstream name;
         name << "z" << i+1 << "_" << j+1;
         std::string name2;
         name >> name2;
         z2[i][j] = modelprice->addVar(0.0,GRB_INFINITY,1,GRB_BINARY,name2);// temp objective coefficient, will be updated in loop
         expp += z2[i][j];
         ++currentconst;
         vbasesp.push_back(0);
      }
      modelprice->addConstr(expp == 1);
      ++totalconstraints;
   }
   std::cout << std::endl;
   //Now repeat for 4-indexed variables
   int constrcol = 0;
   for (int i = 0; i < Pnum-1; ++i)
   {
      for (int j = i+1; j < Pnum; ++j)
      {
         for (int k = 0; k < Psnum[i]; ++k)
         {
            for (int m = 0; m < Psnum[j]; ++m)
            {
               std::stringstream name;
               name << "z" << i+1 <<"_" <<k+1 << "_" <<j+1 << "_" <<m+1;
               std::string name2;
               name >> name2;
               double prod = 2*lambda[i]*lambda[j]*(PsuppT[i][k]*PsuppT[j][m]);
               z4[i][k][j][m] = modelprice->addVar(0.0,GRB_INFINITY,prod,GRB_CONTINUOUS,name2);
               expp2[constrcol]+= z4[i][k][j][m]-z2[i][k];
               ++constrcol;
               expp2[constrcol]+= z4[i][k][j][m]-z2[j][m];
               ++constrcol;
               vbasesp.push_back(0);
            }
         }
      }
   }

   //Add remaining constraints to model
   for (int k = 0; k < constrcol; ++k)
   {
      modelprice->addConstr(expp2[k] <= 0);
      ++totalconstraints;
   }

   //Now solve: IP, then master - should be master then IP?
   double * yhat = model->get(GRB_DoubleAttr_Pi, model->getConstrs(), totalsupp);

   // add/update the coefficient on z[i][j]
   currentconst = 0;
   for (int i = 0; i < Pnum; ++i)
   {
      for (int j = 0; j < Psnum[i]; ++j)
      {
         z2[i][j].set(GRB_DoubleAttr_Obj,yhat[currentconst]-zikcost[i][j]-zjmcost[i][j]);
         ++currentconst;
      }
   }
   modelprice->optimize();
   //To write model to a file
/*   std::ostringstream filename2;
   filename2 << "/Users/spatterson/ForXcode/PricingLP/Pricing_test_inCG" << iter <<".lp";
   modelprice->write(filename2.str());*/
   
   int maxiter = 8;//Pnum*10;
   bool isneg = true;
   while ( isneg )
   {

      if (iter >= maxiter)
      {
         break;
      }
      isneg = false;

      //Create column of correct coefficients for constraints in master problem
      double expcoeff2[totalsupp];
      int l = 0;
      long int h = 0;
      long int Pprod = S0;
      for (int i = 0; i < Pnum; ++i)
      {
         Pprod /= Psnum[i];
         for (int j = 0; j < Psnum[i]; ++j)
         {
            if (z2[i][j].get(GRB_DoubleAttr_X) > 0)
            {
               expcoeff2[l] = 1;
               h += Pprod*j;
            }
            else
            {
               expcoeff2[l] =0 ;
            }
            ++l;
         }
      }
      if (colin[h] == 1)
      {
         std::cout << "Reintroducing same column." << std::endl;
         break;
      }
      
      //Add new column to master problem
      GRBColumn newcol;
      newcol.addTerms(expcoeff2, constrs, totalsupp);
      w.push_back(model->addVar(0.0, GRB_INFINITY, c[h], GRB_CONTINUOUS, newcol));//This updates the objective
      vbases.push_back(-1);
      colin[h] = true;
      if (warmstart == 1 )
      {
         model->update();
         for (int i = 0; i < w.size(); ++i)
         {
            w[i].set(GRB_IntAttr_VBasis, vbases[i]);
         }
         for (int i = 0; i < totalsupp; ++i)
         {
            constrs[i].set(GRB_IntAttr_CBasis, cbases[i]);
         }
      }
      vbstart = vbases;
      oldobj = newobj;
      model->optimize();
//         std::ostringstream fname2;
//         fname2 <<"/Users/spatterson/ForXcode/Barycenter/Column Generation/TestNewGen_"<<iter<<".lp";
//         model->write(fname2.str());
      newobj = model->get(GRB_DoubleAttr_ObjVal);

      ++iter;
      if (warmstart == 1)
      {
         for (int i = 0; i < w.size(); ++i)
         {
            vbases[i] = (w[i].get(GRB_IntAttr_VBasis));
//               std::cout << vbases[i] <<std::endl;
         }
         for (int i = 0; i < totalsupp; ++i)
         {
            cbases[i] = constrs[i].get(GRB_IntAttr_CBasis);
         }
      }
      
      
      t0 = t;
      t = std::chrono::steady_clock::now();

   //      double currentobj = model->get(GRB_DoubleAttr_ObjVal);
   //      std::cout << currentobj << std::endl;
      double * yhat = model->get(GRB_DoubleAttr_Pi, model->getConstrs(), totalsupp);

      // add/update the coefficient on z[i][j]
      currentconst = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         for (int j = 0; j < Psnum[i]; ++j)
         {
            z2[i][j].set(GRB_DoubleAttr_Obj,yhat[currentconst]-zikcost[i][j]-zjmcost[i][j]);
            if (pwarmstart == 1)
            {
               z2[i][j].set(GRB_DoubleAttr_Start, vbasesp[currentconst]);
            }
            ++currentconst;
         }
      }
/*      if (pwarmstart == 1)
      {
         for (int i = 0; i < Pnum-1; ++i)
         {
            for (int j = i+1; j < Pnum; ++j)
            {
               for (int k = 0; k < Psnum[i]; ++k)
               {
                  for (int m = 0; m < Psnum[j]; ++m)
                  {
                     z4[i][k][j][m].set(GRB_DoubleAttr_Start, GRB_UNDEFINED);
                     ++currentconst;
                  }
               }
            }
         }
      }*/
      modelprice->optimize();
      //To write model to a file
/*      std::ostringstream filename2;
      filename2 << "/Users/spatterson/ForXcode/PricingLP/Pricing_test_inCG" << iter <<".lp";
      modelprice->write(filename2.str());*/
      
      t0 = t;
      t = std::chrono::steady_clock::now();
      pricingtimes.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(t-t0).count() );
      if (pwarmstart == 1)
      {
         int tempind = 0;
         for (int i = 0; i < Pnum; ++i)
         {
            for (int j = 0; j < Psnum[i]; ++j)
            {
               vbasesp[tempind] =z2[i][j].get(GRB_DoubleAttr_X);
               ++tempind;
            }
         }
/*         for (int i = 0; i < Pnum-1; ++i)
         {
            for (int j = i+1; j < Pnum; ++j)
            {
               for (int k = 0; k < Psnum[i]; ++k)
               {
                  for (int m = 0; m < Psnum[j]; ++m)
                  {
                     vbasesp[tempind] = (z4[i][k][j][m]).get(GRB_DoubleAttr_X);
                     ++tempind;
                  }
               }
            }
         }*/
/*         for (int i = 0; i < totalconstraints; ++i)
         {
            cbases[i] = constrs[i].get(GRB_IntAttr_CBasis);
         }*/
      }
      
      //Need to keep thinking about correct requirement for continuing.
      //If let maxiter dictate, end up re-introducing the same column.
      if (modelprice->get(GRB_DoubleAttr_ObjVal) >1e-4)
      {
         isneg = true;
         /*std::cout << "Solution to IP:" <<std::endl;
         for (int i = 0; i < Pnum; ++i)
         {
            for (int j = 0; j < Psnum[i]; ++j)
            {
               std::cout << z2[i][j].get(GRB_DoubleAttr_X) <<std::endl;
            }
         }*/
      }
   }

   t = std::chrono::steady_clock::now();
   std::cout << "Total Running time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-start).count()  <<"ms" << std::endl;
   std::cout << "Total Columns Introduced: " << w.size() << std::endl;
   std::cout << "Optimal value: " << newobj << " found after " << iter << " solves." << std::endl;
   std:: cout << "Average Pricing Time " << std::accumulate(pricingtimes.begin(), pricingtimes.end(),0)/iter << "ms" <<std::endl;

//Tracking fractional solutions only makes sense for the LP relaxation
//   std::cout << "Number of fractional LP solutions is " << fracsol << " out of " << iter << " iterations." <<std::endl;
   
   std::cout << "Enter to end." <<std::endl;
   std::cin.get();
   model->terminate();
   model->reset();
   modelprice->terminate();
   delete modelprice;

   delete model;
   delete env;
   return 0;
}

   bool validateDist(const std::vector<SuppPt>::iterator it1, const std::vector<SuppPt>::iterator it2)
   {
      std::vector<SuppPt>::iterator it = it1;
      double total = 0;
      while (it != it2)
      {
         total += it->mass;
         ++it;
      }
      
      if (total < 0.95 || total > 1.05)
      {
         std::cout << "Warning: Far from 1. Consider alternative." << total << std::endl;
      }
      if (total != 1)
      {
         it1->mass = it1->mass + (1.0-total);
      }
      return true;

}

int makegreedy(std::vector<SuppPt> &Psupp2, std::vector<int> &Psnum2, int * Psnum, const int &Pnum, std::vector<double> &winit,const long int &S0)
{
   std::vector<int> indices(Pnum,0);
   std::vector<int> startindices(Pnum,0);
   
   for (int i = 1; i < Pnum; ++i)
   {
      indices[i] = indices[i-1];
      indices[i] += Psnum2[i-1];
      startindices[i] = indices[i];
   }
   
   while (Psupp2[Psnum2[0]-1].mass >1e-12)
   {
      //      std::cout << Psupp2[Psnum2[0]-1].mass <<std::endl;
      double minmass = 1;
      for (int i = 0; i < Pnum; ++i)
      {
         if (Psupp2[indices[i]].mass < minmass)
         {
            minmass = Psupp2[indices[i]].mass;
         }
      }
      
      long int loocur = S0;
      int unki = 0;
      int index = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         unsigned long int index2 = Psupp2[indices[i]].combos-index;
         loocur /= Psnum[i];
         unki += loocur*index2;
         index += Psnum[i];
         Psupp2[indices[i]].mass -= minmass;
         if (Psupp2[indices[i]].mass <= 1e-12)
         {
            ++indices[i];
         }
      }
      
      winit[unki] = minmass;
   }
   return 0;
}


std::vector<double> dualtotals(const double * const yhat, const int Amrn, const long int blocklength, const int * Psnum, const int Pnum)
{
   std::vector<double> ytotals(blocklength,0);
   int i =0; //Should really pass this in. Related to numinprice
   int totalPs = 0;
   long int startindex = 0;
   long int blockinum = 1;
   long int blockilength = blocklength;
   long int blocki1length = blockilength/Psnum[i];
   long int j = 0;
   long int index = 0;
   while (j < Amrn)
   {
      if (yhat[j] != 0)
      {
         //         std::cout << yhat[j] << std::endl;
         index = 0;
         for (int k = 0; k < blockinum; ++k)
         {
            for (int l = 0; l < blocki1length; ++l)
            {
               //               std::cout << "index: " << startindex+index << " yhat: " << yhat[j] << std::endl;
               ytotals[startindex+index] += yhat[j];
               ++index;
            }
            index = blockilength*(k+1);
         }
      }
      ++j;
      //      std::cout << j << std::endl;
      startindex += blocki1length;
      if (j-totalPs == Psnum[i])
      {
         blockinum *= Psnum[i];
         blockilength /= Psnum[i];
         startindex = 0;
         totalPs += Psnum[i];
         ++i;
         if (i < Pnum)
         {
            blocki1length /= Psnum[i];
         }
      }
   }
   return ytotals;
}
