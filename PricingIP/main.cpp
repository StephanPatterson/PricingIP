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
double compute_c(const long int &, const int &, const int *, const std::vector< std::vector<SuppPt> > &, const long int &, const double *);

int main(int argc, const char * argv[]) {
   
   auto start = std::chrono::steady_clock::now();
   auto t = start;
   auto t0 = t;
   std::vector< std::time_t > pricingtimes;

   std::cout.precision(10);
   int Pnum =14;//If NumMoMeas = 1, Number of months to process. Will go consecutively from start month, skipping any empty months.
   std::string startmonth = "Sep";
   int startyear = 14;
   bool namevars = 0;
   bool checkcolfordup = 0;

   int NumMoMeas = 1;//Number of Months in a single Measure. Use 1 for original approach of each month is its own measure
   //Otherwise, Number of measures that will be created. Will combine NumMoMeas of months into each measure
//   int fracsol = 0;
   srand(time(NULL));
   
   int warmstart = 1; //If 1, the previous solution is saved and reintroduced before ->optimize()
   // This should be 1 for fastest running times
   int pwarmstart = 1; //Separate warm-start toggle for modelprice->optimize.
   
   int Psnum[Pnum];
   //For the new IP, it makes more sense to have a two-indexed vector PsuppT instead of one long vector Psupp.
   std::vector< std::vector<SuppPt> > PsuppT(Pnum);
   int totalsupp = 0;
   int indices[Pnum];
   double lambda[Pnum];

   std::string temp;

   std::ifstream indata;
   std::ostringstream fname;
   int year;
   std::string month;
   if (argc > 1)
   {
      fname << argv[1];
      indata.open(fname.str());
   }
   else
   {
      std::cout << "Please enter filepath for DenverCrime file." << std::endl;
      std::string filename;
      getline( std::cin, filename);
      indata.open(filename);
   }

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
         std::cout << "Invalid Input. End of file reached while searching for starting Month/Year." << std::endl;
         return 1;
      }
   }
   
   double loc1;
   double loc2;
   indata >> loc1 >> loc2;
   std::string currentmonth = startmonth;
         
   for (int i = 0; i < Pnum; ++i)
   {
      Psnum[i] = 0;
      indices[i] = 0;
      for (int j = 0; j < NumMoMeas; ++j)
      {
         while (month == currentmonth )
         {
            //adding shift so that all coordinates are positive
            PsuppT[i].push_back(SuppPt(loc1+115, loc2, 1.0,totalsupp));
            ++Psnum[i];
            ++totalsupp;
                  
            indata >> year >> month >> loc1 >> loc2;
            if (indata.eof())
            {
               indata.close();
               std::cout << "Invalid Input. End of file reached while reading in months to measures." << std::endl;
               return 1;
            }
         }
         currentmonth = month;
      }
      std::cout << "Month measure: " << i+1 << " Size: " << Psnum[i] << std::endl;
      
      //Scale the masses for each month to sum to 1
      double totalmass = 0;
      for (int j = 0; j < Psnum[i]-1; ++j)
      {
         PsuppT[i][j].mass /= Psnum[i];
         totalmass += PsuppT[i][j].mass;
      }
      PsuppT[i][Psnum[i]-1].mass = 1-totalmass;
      currentmonth = month;
   }
   indata.close();
   
   //Compute the weights for each month
   double sum = 0;
   for (int i = 0; i < Pnum-1; ++i)
   {
      lambda[i] = (double)Psnum[i]/totalsupp;
      sum += lambda[i];
   }
   lambda[Pnum-1] = 1-sum;

   //Calculate number of possible supp pts S0
   long int S0 = 1;
   for (int i = 0; i < Pnum; ++i)
   {
      S0 *= Psnum[i];
   }
   std::cout << "Size of S0: " << S0 << std::endl;
   
   std::vector<bool> colin;
   if (checkcolfordup)
   {
      //For ensuring a column does not get added repeatedly
      colin.resize(S0,false);
   }

   //Model for Master Problem
   GRBEnv* env = new GRBEnv();
   GRBModel* model = new GRBModel(*env);
   model->set(GRB_IntParam_Method, 0);
   model->set(GRB_IntParam_Presolve, 0);
   model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
   model->set(GRB_IntParam_OutputFlag, 0); //Turning off display on screen. Disable to see iterations, objective value, etc.
   
   std::vector< GRBVar > w;
   GRBLinExpr  exp[totalsupp];
   std::vector<int> vbases;
   double oldobj;
   double newobj;
   
   //Create a temporary copy of the set of support points, will be modified
   std::vector< std::vector<SuppPt> >::iterator Psuppit = PsuppT.begin();
   std::vector<std::vector<SuppPt> > Psupp2(Pnum);
   std::copy(Psuppit, Psuppit+Pnum, Psupp2.begin() );
   //Using index as a dummy variable for readability
   int index = 0;
   
   while (Psupp2[0][Psnum[0]-1].mass >1e-15)
   {
      //Find smallest remaining mass among support points of current combination
      double minmass = 1;
      for (int i = 0; i < Pnum; ++i)
      {
         if (Psupp2[i][indices[i]].mass < minmass)
         {
            minmass = Psupp2[i][indices[i]].mass;
         }
      }
//      std::cout << "Minimum mass is " << minmass <<std::endl;
      
      //Index math based on Algorithm 1 in A Column Generation Approach to the Discrete Barycenter Problem
      long int denom = S0;
      int h = 0;
      int index = 0;
      
      double sum1 = 0;
      double sum2 = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         sum1 += lambda[i]*PsuppT[i][indices[i]].loc1;
         sum2 += lambda[i]*PsuppT[i][indices[i]].loc2;

         denom /= Psnum[i];
         unsigned long int index2 = Psupp2[i][indices[i]].combos-index;
         h += denom*index2; //computing h for later computation of cost c
         index += Psnum[i];
         Psupp2[i][indices[i]].mass -= minmass;
      }
      
      //Create a variable for this combination
      double newc = 0;
      SuppPt Pbar0 = SuppPt(sum1, sum2, 0.0);
      for (int i = 0; i < Pnum; ++i)
      {
         double diff1 =Pbar0.loc1-PsuppT[i][indices[i]].loc1;
         double diff2 =Pbar0.loc2-PsuppT[i][indices[i]].loc2;
         newc += lambda[i]*(diff1*diff1+diff2*diff2);
         //Updating indices needed to be done later than previous code projects to allow for cost to be computed here
         if (Psupp2[i][indices[i]].mass <= 1e-15)
         {
            ++indices[i];
         }
      }
      
      w.push_back(model->addVar(0.0, GRB_INFINITY, newc, GRB_CONTINUOUS));
      vbases.push_back(0);
      
      //Add variable to corresponding constraints
      denom = S0/Psnum[0];
      int startindex = Psnum[0];
      int jindex = floor(h/denom);
      exp[jindex] += *(w.end()-1);
      for (int l = 1; l < Pnum-1; ++l)
      {
         denom /= Psnum[l];
         jindex = startindex+floor( (h % (denom*Psnum[l]))/denom);
         exp[jindex] += *(w.end()-1);
         startindex += Psnum[l];
      }
      jindex = startindex+ h % Psnum[Pnum-1];
      exp[jindex] += *(w.end()-1);
   }

   index = 0;
   for (int i = 0; i < Pnum; ++i)
   {
      for (int j = 0; j < Psnum[i]; ++j)
      {
         model->addConstr(exp[index] == PsuppT[i][j].mass);
         ++index;
      }
   }
   
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
   std::vector< std::vector< GRBVar >> z2(Pnum);
   std::vector< std::vector< std::vector< std::vector< GRBVar >> >>z4(Pnum,std::vector< std::vector< std::vector< GRBVar >>>(Psmax, std::vector< std::vector< GRBVar >>(Pnum)));
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
   //Add two-indexed variables to the IP
   int currentconst = 0;
   int totalconstraints = 0;
   for (int i = 0; i < Pnum; ++i)
   {
      GRBLinExpr expp;
      for (int j = 0; j < Psnum[i]; ++j)
      {
         if (namevars)
         {
            std::stringstream name;
            name << "z" << i+1 << "_" << j+1;
            std::string name2;
            name >> name2;
            z2[i].push_back( modelprice->addVar(0.0,GRB_INFINITY,1,GRB_BINARY,name2));// temp objective coefficient, will be updated in loop
         }
         else
         {
            z2[i].push_back( modelprice->addVar(0.0,GRB_INFINITY,1,GRB_BINARY));
         }
         expp += z2[i][j];
         ++currentconst;
         vbasesp.push_back(0);
      }
      modelprice->addConstr(expp == 1);
      ++totalconstraints;
   }
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
               if (namevars)
               {
                  std::stringstream name;
                  name << "z" << i+1 <<"_" <<k+1 << "_" <<j+1 << "_" <<m+1;
                  std::string name2;
                  name >> name2;
                  double prod = 2*lambda[i]*lambda[j]*(PsuppT[i][k]*PsuppT[j][m]);
                  z4[i][k][j].push_back( modelprice->addVar(0.0,GRB_INFINITY,prod,GRB_CONTINUOUS,name2));
                  
               }
               else
               {
                  double prod = 2*lambda[i]*lambda[j]*(PsuppT[i][k]*PsuppT[j][m]);
                  z4[i][k][j].push_back( modelprice->addVar(0.0,GRB_INFINITY,prod,GRB_CONTINUOUS));
               }
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
   
   int maxiter = 10; //Pnum*5;
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
      if (checkcolfordup)
      {
         if (colin[h] == 1)
         {
            std::cout << "Reintroducing same column." << std::endl;
            break;
         }
      }
      
      //Add new column to master problem
      GRBColumn newcol;
      newcol.addTerms(expcoeff2, constrs, totalsupp);
      
      //Compute newc here
      double newc = compute_c(S0, Pnum, Psnum, PsuppT, h, lambda);
      w.push_back(model->addVar(0.0, GRB_INFINITY, newc, GRB_CONTINUOUS, newcol));//This updates the objective
      vbases.push_back(-1);
      if (checkcolfordup)
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

      //This option reduces the amount of memory used by Gurobi, but slows solution times significantly
      modelprice->reset(1);
      
      
      
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

      /* //No discernable change in solving speeds with this reset
      if (pwarmstart == 1)
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
      
      long int denom = S0;
      int h = 0;
      int index = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         unsigned long int index2 = Psupp2[indices[i]].combos-index;
         denom /= Psnum[i];
         h += denom*index2;
         index += Psnum[i];
         Psupp2[indices[i]].mass -= minmass;
         if (Psupp2[indices[i]].mass <= 1e-12)
         {
            ++indices[i];
         }
      }
      
      winit[h] = minmass;
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

//Function for computing the cost c associated with a particular combination
//Implementation of Algorithm 1 from Column Generation paper
//Inputs:
//S0: product of sizes of measures
//Pnum: number of measures
//Psnum: number of support points per measure
//Psupp: vector containing all support points
//indices: vector containing indices of
double compute_c(const long int &S0, const int &Pnum, const int * Psnum, const std::vector< std::vector<SuppPt> > &Psupp, const long int &h, const double * lambda)
{
   int index = 0;
   int indices[Pnum];
   long int denom = S0/Psnum[0];
   index = floor(h/denom);
   indices[0] = index;

   double sum1 = 0;
   double sum2 = 0;
   sum1 += lambda[0]*Psupp[0][index].loc1;
   sum2 += lambda[0]*Psupp[0][index].loc2;
   for (int i = 1; i < Pnum-1; ++i)
   {
      denom /= Psnum[i];
      index = floor((h%(denom*Psnum[i]))/denom);
      indices[i] = index;

            
      sum1 += lambda[i]*Psupp[i][index].loc1;
      sum2 += lambda[i]*Psupp[i][index].loc2;
   }
   index = h%Psnum[Pnum-1];
   indices[Pnum-1] = index;
   sum1 += lambda[Pnum-1]*Psupp[Pnum-1][index].loc1;
   sum2 += lambda[Pnum-1]*Psupp[Pnum-1][index].loc2;

   //Create a variable for this combination
   double newc = 0;
   SuppPt Pbar0 = SuppPt(sum1, sum2, 0.0);
   for (int i = 0; i < Pnum; ++i)
   {
      double diff1 =Pbar0.loc1-Psupp[i][indices[i]].loc1;
      double diff2 =Pbar0.loc2-Psupp[i][indices[i]].loc2;
      newc += lambda[i]*(diff1*diff1+diff2*diff2);
   }
   return newc;
}
