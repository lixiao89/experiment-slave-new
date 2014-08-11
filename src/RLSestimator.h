#ifndef _RLS_ESTIMATOR_H
#define _RLS_ESTIMATOR_H

#include<iostream>
#include<cmath>

class RLSestimator{

    public:

        // x[0] = mu is the cofficient of kinetic friction, x[1] = Fc is the estimated cutting force
        vctFixedSizeVector<double,2> x;  
        vctFixedSizeVector<double,2> xlast; // used for motionDetection feature
        // P is the covariance of x        
        vctFixedSizeMatrix<double,2,2> P;
        // Fest is the estimated total tangential force 
        double Fest;
        double FestLast; // used for motionDetection feature

        // indicate a cutting failure mode
        bool fail;
        // indicate contate state
        bool contactState;
        // true if wam is not moving;
        bool wamMotionState;
        
        // SC: sliding and cutting
        // SWC: sliding without cutting
        // CWS: cutting without sliding
        enum scenerio {IDLE, SC, SWC, CWS};
        scenerio sc;

        // error threshold
        double cuttingFailureThreshold; // defines when cutting failure happends
        double contactThreshold;//

        //system matrix yk = Hk.transpose*x + v
        vctFixedSizeMatrix<double,2,2> Hk;
        // measurement error covariance
        vctFixedSizeMatrix<double,2,2> Rk;
        
        // yk is the measurement vector with first element measured tangential force and second element 0
        vctFixedSizeVector<double,2> yk;
        // constructor
        RLSestimator(vctFixedSizeVector<double,2>& xinit): 
        
         x(xinit),
         xlast(xinit),
         P(vct2x2::Eye()),
         Rk(vct2x2::Eye()),
         fail(false),
         contactState(true),
         wamMotionState(true), 
         sc(SWC),  
         cuttingFailureThreshold(5),
         contactThreshold(0.01),
         Hk(vct2x2::Eye()),
         Fest(0),
         FestLast(0),
         yk(0,0){
                
               Hk[0][0] = 0;
               Hk[0][1] = 0;
               Hk[1][0] = 1;
               Hk[1][1] = 0;
         }

        ~RLSestimator(){};
        
        //returns true if the cutter is in contact with the cutting surface
        bool inContact(const double& Fn){
            if(std::abs(Fn) > contactThreshold && contactState == true){
                std::cout<< "Cutter in Contact" << std::endl;
                contactState = false;
                return true;
            }
            if(std::abs(Fn) <= contactThreshold && contactState == false){
                std::cout<< "Cutter not in Contact!" << std::endl;
                contactState = true;
                return false;
            }
        }

        // detects if WAM is in motion
        bool WAMnotMoving(const bool& wamNotMoving){
            if(wamNotMoving && wamMotionState == true){
                std::cout<< "WAM STOPPED!" << std::endl;
                wamMotionState = false;
            }
            if(!wamNotMoving && wamMotionState == false){
                std::cout<< "WAM moving" << std::endl;
                wamMotionState = true;
            }
        }       

        void GetEstimates (vctFixedSizeVector<double,2>& xesti, double& Festi) const{

            xesti[0] = x[0];
            xesti[1] = x[1];
            Festi = Fest; 
        }        
        
        // returns true if detects a cutting failure
        // parameters:
        //
        // Fn: measured normal force
        // Fe: measured tangential force
        // P: covariance matrix from last step 
        // xnew: newly estimated [mu, Fc]
        // Fnew: newly estimated tangential force
        // Cov: updated covariance matrix
        // isCuting: true if sliding and cutting, false if sliding but not cutting
       void Evaluate( const double &Fn, 
                      const double &Fe, 
                      const vctFixedSizeMatrix<double,2,2>& p, 
                      vctFixedSizeVector<double,2>& xnew, 
                      double& Fnew, 
                      vctFixedSizeMatrix<double,2,2>& Cov, 
                      bool isCutting){

               // xnew = xlast;
               // Fnew = Flast;
                
                // in the case of sliding without cutting
                if(!isCutting){
                    Hk[1][0] = 0;
                    xnew[1] = 0;
                }

                Hk[0][0] = Fn;
                yk[0] = Fe;
                
                vctFixedSizeMatrix<double,2,2> K;
                vctFixedSizeMatrix<double,2,2> tempK;
                vctFixedSizeMatrix<double,2,2> invtempK;

                tempK = Hk.Transpose()*P*Hk + Rk;

                Inverse(tempK,invtempK);
                K = P*Hk*invtempK;

                xnew = xnew + K*(yk - Hk.Transpose()*x);

                Cov = (vct2x2::Eye() - K*Hk.Transpose())*P;
                
                Fnew = xnew[0]*Fn + xnew[1];
             
                //xlast = xnew;
                //Flast = Fnew;

               }


        // Evaluate the estimator comparing the outcome of three different cases
        // 1. sliding and cutting (simulataneously estimating mu and Fc)
        // 2. sliding without cutting (Fc = 0)
        bool EvaluateWithComparison( const double &Fn, const double &Fe ){

                vctFixedSizeVector<double,2> x1(xlast);
                vctFixedSizeVector<double,2> x2(xlast);

                double F1 = FestLast;
                double F2 = FestLast;

                vctFixedSizeMatrix<double,2,2> Cov1;
                vctFixedSizeMatrix<double,2,2> Cov2;

                // evaluate the first scenerio
                Evaluate(Fn, Fe, P, x1, F1, Cov1,true);
                // evaluate the second scenerio
                Evaluate(Fn, Fe, P, x2, F2, Cov2, false);
                
                double Fdiff[2] = {std::abs(F1-Fe), std::abs(F2-Fe)};

                // finding the case that best resembles the measured tangential force
                int minIndex = 0;
                double minTemp = 100;
                for(int j = 0; j < 2; ++j){

                    if(Fdiff[j] < minTemp){
                        minIndex = j;
                        minTemp = Fdiff[j];
                    }
                }
                
                // find the most suitable scenerio and update member variables
                if(minIndex == 0){
                    x = x1;
                    Fest = F1;
                    P = Cov1;
                }
                else{
                    x = x2;
                    Fest = F2;
                    P = Cov2;
                }

                // for printing out wam status
                if(minIndex == 0 && sc == SC){
                    std::cout<< "Cutting and sliding" << std::endl;
                    sc = SWC;
                }
                if(minIndex == 1 && sc == SWC){
                    std::cout<< "sliding without cutting" << std::endl;
                    sc = SC;
                }
                
                xlast = x;
                FestLast = Fest;

              if(fabs(Fest - Fe) > cuttingFailureThreshold){

                return true;
            }
            else{

                return false;
            }


        }





       /** parameters: 
        *
        * Fn: normal force as measured by the force sensor (currently force in the z direction)
        * Fe: tangential force as measured by the force sensor (currently force in the x direction)
        * currtime: current time = osaGetTime() - startTime
        * currJointPos: current joint positions for all 7 joints
        * rlsEstData: is a vector containing [Fe, Fn, mu, Fc, Fest, haveFailed] that gets pushed to ROS (gets changed by the method)
        * failstate: used to indicate cutting failure (true if failed)  (gets changed by the method)
        * haveFailed: a double value used to draw a spike on rqt if failure happens (set to a large number)  (gets changed by the method)
        * wamNotMoving: is true when wam is not moving
        * motionDetection: a switch to toggle between with and without motion detection feature
        */

        
   void RLSestimate(const double &Fn, 
                    const double &Fe, 
                    const double& currtime,
                    const vctDynamicVector<double>& currJointPos, 
                    vctDoubleVec& rlsEstData,
                    bool& failstate,
                    double& haveFailed,
                    const bool& wamNotMoving,
                    const std::string motionDetection = "OFF"){

        // detect is cutter is in contact
        inContact(Fn);
        // detects if wam is motionDetection
        WAMnotMoving(wamNotMoving);

        if (motionDetection == "OFF"){

                        if(this->EvaluateWithComparison( Fn, Fe ) && !failstate){
                   std::cout<<"Cutting Failure at time:"<< currtime <<std::endl;
                   failstate = true;
                   haveFailed = 70;
                }
                else{
                    failstate = false;
                    haveFailed = 0;
                }

                rlsEstData.Assign((double)Fe, Fn, x[0], x[1], Fest, haveFailed);
        }
        if (motionDetection == "ON"){

                
            if(wamNotMoving){

                x.Assign((double)0,0);
                Fest = 0;
                haveFailed = 0;
                rlsEstData.Assign((double)Fe, Fn, x[0], x[1], Fest, haveFailed);
      
               }
              else{
                if(this->EvaluateWithComparison( Fn, Fe ) && !failstate){
                   std::cout<<"Cutting Failure at time:"<< currtime <<std::endl;
                   failstate = true;
                   haveFailed = 70;
                }
                else{
                    failstate = false;
                    haveFailed = 0;
                }

                rlsEstData.Assign((double)Fe, Fn, x[0], x[1], Fest, haveFailed);
            }
    }
 }

      void Inverse(vctFixedSizeMatrix<double,2,2>& M, vctFixedSizeMatrix<double,2,2>& Minv){
            double detM;
            detM = M[0][0]*M[1][1] - M[0][1]*M[1][0];
            
            if(abs(detM) < 0.0001){
                detM = 0.0001;
            }

            Minv[0][0] = M[1][1];
            Minv[0][1] = -M[0][1];
            Minv[1][0] = -M[1][0];
            Minv[1][1] = M[0][0];

            Minv = (1/detM)*Minv;


        }
};




#endif
