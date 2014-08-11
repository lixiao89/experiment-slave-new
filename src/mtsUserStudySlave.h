/*-*- Mode++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-   */
/*ex: set filetype=cpp softtabstop=4 shiftwidth=4 tabstop=4 cindent expandtab:*/

/*
  $Id: mtsUserStudySlave.h 3181 2011-11-15 15:41:28Z sleonar7 $

  Author(s):  Simon Leonard
  Created on: 2013-07-22

  (C) Copyright 2012 Johns Hopkins University (JHU), All Rights Reserved.

--- begin cisst license - do not edit ---

This software is provided "as is" under an open source license, with
no warranty.  The complete license can be found in license.txt and
http://www.cisst.org/cisst/license.txt.

--- end cisst license ---
*/

#include <fstream>
#include <cisstCommon/cmnPath.h>
#include <cisstOSAbstraction/osaSleep.h>
#include <cisstOSAbstraction/osaGetTime.h>

#include <cisstMultiTask/mtsTaskPeriodic.h>
#include <cisstMultiTask/mtsInterfaceRequired.h>
#include <cisstMultiTask/mtsInterfaceProvided.h>

#include <cisstVector.h>
#include <cisstParameterTypes.h>

#include <sawJR3ForceSensor/osaJR3ForceSensor.h>
#include <sawControllers/osaGravityCompensation.h>
#include "osaHybridForcePosition.h"

#include <cisstRobot/robLinearSE3.h>

#include <cisstNumerical/nmrLSMinNorm.h>
#include <cisstNumerical/nmrSavitzkyGolay.h>

#include"RLSestimator.h"
class mtsUserStudySlave : public mtsTaskPeriodic {

private:

  // robot stuff
  robManipulator robot;         // robot kinematics
  robManipulator tool;          // tool kinematics
  robLinearSE3* se3traj;        // interpolation 
  std::list<double> dt;
  double timer;

  vctVec    qsold;
  vctFrm4x4 Rtwtsold; // previous command to the slave
  vctRot3   Rts;      // rotation of sensor wrt tool

  //vctMatrixRotation3<double> Rts;
  
  // JR3 sensor
  osaJR3ForceSensor jr3;

  // the controllers
  osaHybridForcePosition hfp;
  
  // interfaces
  mtsInterfaceProvided* mastertel;  // telemetry (slave->master)
  mtsInterfaceProvided* mastercmd;  // commands  (master->slave)
  mtsInterfaceRequired* slave;      // interface to the slave PIDs
  mtsInterfaceRequired* control;    // control interface

  enum State{ IDLE, ENABLE };
  State state;
  bool IsEnabled() const { return state == ENABLE; }


  // Slave side
  //! Read the joint positions
  mtsFunctionRead  mtsGetPosition;
  //! Write the joint positions
  mtsFunctionWrite mtsSetPosition;

  // Master side 
  prmPositionCartesianGet prmCommandSE3;    // Cartesian command (master side)
  prmForceCartesianGet    prmCommandWrench; // Wrench command (master side)
  prmPositionCartesianGet prmTelemetrySE3;  // Cartesian telemetry
  prmPositionJointGet     prmTelemetryRn;   // Joint telemetry


  //! Get the position command increment from the master
  bool GetCommand( vctFrm4x4& Rtwts );

  //! Get the wrench command from the master
  bool GetCommand( osaJR3ForceSensor::Wrench& ws );

  //! Get the measured joint angles
  bool GetPosition( vctVec& q );

  //! Get the measured cartesian position
  bool GetPosition( vctFrm4x4& Rtwt );

  //! Get the measured wrench
  bool GetWrench( osaJR3ForceSensor::Wrench& w );

  //! Get the position command increment from the master
  void GetCommand();

  //! Set the telemetry 
  void SetTelemetry();

  void IdleMotion();
  void HybridMotion();
  
  vctDynamicVector<double> sg;
  std::list< osaJR3ForceSensor::Wrench > stdft;    

//--------------- For RLS --------------------
  RLSestimator* rls;
  std::ofstream ofsForceData;
  double startTime;
  bool failstate;
  bool isMoving;
  double prevTime;


  vctFrame4x4<double> prevPos;// record the last position manipulator is in
  vctDynamicVector<double> prevJointPos;

  std::vector< std::vector<double> > jointPoses;
  std::vector<double> timeStamps;
  int avgNum;// number of joint position measurement to take in calculation of velocity

  // Estimated quantities from RLSestimator, xesti = [mu, Fc], Festi is the estimated tangential force
  vctFixedSizeVector<double,2> xesti;
  double Festi;
  vctDoubleVec rlsEstData;

  mtsInterfaceProvided* rlsProvided; // used to connect to mtsROSBridge;
  double haveFailed; // for cisst_msgs::rlsEstData
 //--------------------------------------------

 public:
  int logcnt;
  double logtime;

  void Toggle(){ 
    if( state == ENABLE ) {
      ofs.close();
      state = IDLE; 
      ofs2.open( "log.txt" );
    }
    else{
      char fname[32];
      sprintf( fname, "log%d.txt", logcnt++ );
      logtime=osaGetTime();
      ofs.open( fname );
      ofs2.close();
      state = ENABLE; 
    }
  }

 public:
  std::ofstream ofs;
  std::ofstream ofs2;
  //!
  /**
     @param name The component name
     @param period The component period (0.01s)
     @param robotfilename The kinematics file of the robot
     @param Rtw0 Orientation/position of the base wrt world
     @param Rtnt Orientation/position of tool wrt robot last link
     @param qinit Initial joints position of the robot robot
     @param qready Ready joints position of the robot
     @param jr3 JR3 F/T sensor
     @param gc Gravity compensation controller
     @param hfp Hybrid FT/Position controller
   */
  mtsUserStudySlave(const std::string& name,
		    double period,
                    
		    const std::string& robotfilename, 
		    const vctFrm4x4& Rtw0 = 
		    vctFrm4x4( vctRot3(  0.0,  0.0, -1.0,
					 0.0,  1.0,  0.0,
					 1.0,  0.0,  0.0,
					 VCT_NORMALIZE ),
			       vct3( 0.0 ) ),
		    
		    const vctFrm4x4& Rtnt =
		    vctFrm4x4( vctRot3(  0.9511,  0.0000,  -0.3090,
					 0.0000,  1.0000,   0.0000,
					 0.3090,  0.0000,   0.9511,
					 VCT_NORMALIZE ),
			       vct3( 0.0, 0.0, 0.15 ) ),
		    
		    const vctVec& qinit = 
		    vctVec(7.0, 0.0, -cmnPI_2, 0.0, cmnPI, 0.0, -cmnPI_2, 0.0));
    

// --------------------- For RLS -----------------------------------
  void PrintTime(){
      std::cout<< "current time is: "<< osaGetTime() - startTime <<std::endl;
  }


  void CalcVectorAverage(const std::vector<double>& vect, double& avg){

      double temp = 0;
      for(int i = 0; i < vect.size(); i++){
            temp = temp + vect.at(i);
      }

      avg = temp/vect.size();
  }

void CalcAverageVelocity(vctDynamicVector<double>& currJointPos, double& currTime, vctDynamicVector<double>& avgVel){

    std::vector<double> temp;

    for(int i = 0; i < 7; i++){
        temp.push_back(currJointPos[i]);
    }

         jointPoses.push_back(temp);
         timeStamps.push_back(currTime);
        // process only the most recent avgNum data
        if(timeStamps.size() > avgNum){
            jointPoses.erase(jointPoses.begin());
            timeStamps.erase(timeStamps.begin());
        }


        vctDynamicVector<double> currPos( 7 , jointPoses.at(avgNum-1).at(0),jointPoses.at(avgNum-1).at(1),jointPoses.at(avgNum-1).at(2),jointPoses.at(avgNum-1).at(3),jointPoses.at(avgNum-1).at(4),jointPoses.at(avgNum-1).at(5),jointPoses.at(avgNum-1).at(6));

        vctDynamicVector<double> pastPos( 7 , jointPoses.at(0).at(0),jointPoses.at(0).at(1),jointPoses.at(0).at(2),jointPoses.at(0).at(3),jointPoses.at(0).at(4),jointPoses.at(0).at(5),jointPoses.at(0).at(6));


        vctDynamicVector<double> jointPosdiff = (currPos - pastPos).Abs();
        
        double timediff = (timeStamps.at(avgNum-1) - timeStamps.at(0));
        avgVel = jointPosdiff / timediff;

}

bool WAMIsNotMoving( vctDynamicVector<double>& currJointPos, double& currTime){
      
    vctDynamicVector<double> jointAvgVel;
    CalcAverageVelocity(currJointPos, currTime, jointAvgVel);
      

            double motionThreshold = 0.028;

        if(jointAvgVel[0] < motionThreshold && jointAvgVel[1] < motionThreshold && jointAvgVel[2] < motionThreshold && jointAvgVel[3] < motionThreshold && jointAvgVel[4] < motionThreshold && jointAvgVel[5] < motionThreshold && jointAvgVel[6] < motionThreshold){

            //std::cout<<jointAvgVel<<std::endl;
            return true;
        }
        else{

            //std::cout<<jointAvgVel<<std::endl;
            return false;
        }
}
  //----------------------------------------------------------------
 



  void Configure( const std::string& );
  void Startup();
  void Run();
  
};

