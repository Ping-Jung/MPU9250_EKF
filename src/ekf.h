#ifndef EKF_H
#define EKF_H
#endif

#include "Arduino.h"
#include "math.h"
#include "tinyutil.h"

typedef struct MAT{
    FTYPE q[4];                // quaternion 
    FTYPE w[3];                // angular velocity
    FTYPE W[12];             // for calculating the noise of gyroscope
    FTYPE F[16];             // the Jacobian of f
    FTYPE P[16];             // covariance Matrix
    FTYPE Q[16];             // process noise covar mat
    FTYPE z[6];                // measurement
    FTYPE h[6];                // measurement model 
    FTYPE v[6];                // v=z-h
    FTYPE H[24];             // the jacobian of h
    FTYPE R[36];             // the noise of measurement
    FTYPE S[36];             //HPH^T+R
    FTYPE S_inv[36];         // S^-1
    FTYPE K[16];             //kalman gain
    FTYPE Sigma[9];          // gryscope bias
    FTYPE duration;          // time duration

   // for temporary computation
    FTYPE temp[40];
    FTYPE temp2[40];
    FTYPE temp3[40];  
}mat;


class ekf
{
    
private:
    struct MAT param;
    FTYPE bax,bay,baz,bmx,bmy,bmz;
    FTYPE wx,wy,wz,dt;
   
 
public:
 /*reference value*/
 FTYPE yaw,pitch,roll;
 int F_ref[16]={0,-1,-2,-3,1,0,3,-2,2,-3,0,1,3,2,-1,0};;
 int W_ref[12]={-1,-2,-3,0,-3,2,3,0,-1,-2,1,0};
 FTYPE I_4by4[16]={1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
 int g[3]={0,0,-1};
 FTYPE mag[3]={0.7943,0,0.6077};         // magnetic inclination angle [costheta 0 sintheta]
   
 void ekf_init(FTYPE bax, FTYPE bay, FTYPE baz, FTYPE bmx, FTYPE bmy, FTYPE bmz){
  // modi note: yaw angle still require improvement
  param={{0}};
  FTYPE declination=0;  // pending, since it is the effect that should consider after from qtn to euler        
  
  roll=atan(bay/baz);
  pitch=atan(bax/sqrt(bay*bay+baz*baz)); 
  
  Serial.print(" roll: ");
  Serial.print(roll); 
  Serial.print(" pitch: ");
  Serial.print(pitch); 


  FTYPE emx=bmx*cos(pitch)+bmz*sin(pitch);             //x
  FTYPE emy=bmx*sin(roll)*sin(pitch)+bmy*cos(roll)-bmz*sin(roll)*cos(pitch);   //y
  yaw=atan2(emy,emx)+declination;

  Serial.print(" yaw: ");
  Serial.println(yaw);  

  Serial.println("the initial attitude..............................");
  // transform the euler angle to quaternion 
  param.q[0]=cos(roll/2)*cos(pitch/2)*cos(yaw/2)+sin(roll/2)*sin(pitch/2)*sin(yaw/2);
  param.q[1]=sin(roll/2)*cos(pitch/2)*cos(yaw/2)-cos(roll/2)*sin(pitch/2)*sin(yaw/2);
  param.q[2]=cos(roll/2)*sin(pitch/2)*cos(yaw/2)+sin(roll/2)*cos(pitch/2)*sin(yaw/2);
  param.q[3]=cos(roll/2)*cos(pitch/2)*sin(yaw/2)-sin(roll/2)*sin(pitch/2)*cos(yaw/2);
  
  ekf_qtn_to_Euler();
  
  // set P as 4x4 identity matrix
  for(int i=0;i<4;i++){
   param.P[i*(5)]=1;
  }  
  
  }   // use accelerometer and magnetometer to get initial quaternion (v)

 /*prediction step*/
 void ekf_prediction(FTYPE wx, FTYPE wy,FTYPE wz, FTYPE dt){
  /*get angular velocity   check unit or resolution*/ 
   param.w[0]=wx;
   param.w[1]=wy;
   param.w[2]=wz;
   param.duration=dt;

   pred_state_update();
   set_F();
   pred_covariance();

 }

/* correction step*/
 void ekf_correction(FTYPE *Acc, FTYPE *Mag){
   /*get norm acceleration and magnetic field check unit and resolution*/  
  norm(Acc);
  norm(Mag);

  param.z[0]=Acc[0];
  param.z[1]=Acc[1];
  param.z[2]=Acc[2];
  param.z[3]=Mag[0];
  param.z[4]=Mag[1];
  param.z[5]=Mag[2];

  set_h_H();
  correction();

 }  // combine measure covariance,kalman gain, state update, covariance update

/*get measurement noise*/
 void ekf_measumement_noise(FTYPE mbax, FTYPE mbay, FTYPE mbaz, FTYPE mbmx, FTYPE mbmy, FTYPE mbmz ){
   param.R[0]=mbax;
   param.R[7]=mbay;
   param.R[14]=mbaz;
   param.R[21]=mbmx;
   param.R[28]=mbmy;
   param.R[35]=mbmz;

 }  // assume each axis is independent 
 void ekf_gyro_noise(FTYPE bx, FTYPE by, FTYPE bz){
  param.Sigma[0]=bx;
  param.Sigma[4]=by;
  param.Sigma[8]=bz; 
  
 } 

/*others*/
 void ekf_qtn_to_Euler(void){
   /*warning: assume quaternion is normalized */ 
    roll=atan2(2*(param.q[1]*param.q[0]+param.q[2]*param.q[3]),1-2*(param.q[1]*param.q[1]+param.q[2]*param.q[2]))*180.0f/PI;
    pitch=asin(2*(param.q[0]*param.q[2]-param.q[1]*param.q[3]))*180.0f/PI;
    yaw=atan2(2*(param.q[0]*param.q[3]+param.q[1]*param.q[2]),1-2*(param.q[2]*param.q[2]+param.q[3]*param.q[3]))*180.0f/PI;
    if (yaw<0){
       yaw+=360;
    }

 } 

/*prediction step*/
 void pred_state_update(void){
  FTYPE bufw,bufx,bufy,bufz;
  FTYPE halft=0.5*param.duration;
  bufw=param.q[0];
  bufx=param.q[1];
  bufy=param.q[2];
  bufz=param.q[3];

  param.q[0]=bufw-halft*(param.w[0]*bufx+param.w[1]*bufy+param.w[2]*bufz);
  param.q[1]=bufx+halft*(param.w[0]*bufw+param.w[2]*bufy-param.w[1]*bufz);
  param.q[2]=bufy+halft*(param.w[1]*bufw-param.w[2]*bufx+param.w[0]*bufz);
  param.q[3]=bufz+halft*(param.w[2]*bufw+param.w[1]*bufx-param.w[0]*bufy);

  FTYPE buffer=InvSqrt(param.q[0]*param.q[0]+param.q[1]*param.q[1]+param.q[2]*param.q[2]+param.q[3]*param.q[3]);
  param.q[0]=param.q[0]*buffer;
  param.q[1]=param.q[1]*buffer;
  param.q[2]=param.q[2]*buffer;
  param.q[3]=param.q[3]*buffer;
  
  for(int i=0;i<12;i++){
    int temp=abs(W_ref[i]);
    if(temp==1){
      param.W[i]=halft*(param.q[1]);
    }
    else if (temp==2){
      param.W[i]=halft*(param.q[2]);
    }
    else if (temp==3){
       param.W[i]=halft*(param.q[3]);
    }
    else
      param.W[i]=halft*param.q[0];

    if(W_ref[i]<0){
      param.W[i]=param.W[i]*-1;
    }  
  }

 } 
 void set_F(void){
   FTYPE halft=0.5*param.duration;
   for(int i=0;i<16;i++){
    int temp=abs(F_ref[i]);
    if(temp==1){
      param.F[i]=halft*(param.w[0]);
    }
    else if (temp==2){
      param.F[i]=halft*(param.w[1]);
    }
    else if (temp==3){
       param.F[i]=halft*(param.w[2]);
    }
    else
      param.F[i]=1;

    if(F_ref[i]<0){
      param.F[i]=param.F[i]*-1;
    }  
   } 
 }  
 void pred_covariance(void){
  mat_trans(param.W,param.temp2,4,3);   //W^T
  mat_multiply(param.W,param.Sigma,param.temp3,4,3,3);  //W*(Sigma)
  mat_multiply(param.temp3,param.temp2,param.Q,4,3,4);   //Q=W*(Sigma)*(W^T)
  //Serial.println("Q");
  //print_mat(param.Q,4,4);

  mat_trans(param.F,param.temp2,4,4);  //F^T
  mat_multiply(param.F,param.P,param.temp3,4,4,4);    //F*P
  mat_multiply(param.temp3,param.temp2,param.P,4,4,4);  //F*P*F^T
  //Serial.println("FPF^t not add Q");
  //print_mat(param.P,4,4);


  //mat_add(param.Q, param.temp, param.P,16);  //F*P*F^T+Q


 }

/* correction step*/
 void set_h_H(void){
  // set H and h for accelerometer and magnetometer respectively
  for(int j=0;j<12;j++){
    switch (j)
    {
    case 0:
     param.h[j]=2*(param.q[1]*param.q[3]-param.q[0]*param.q[2]);
     param.H[j]=-2*param.q[2];
     param.h[j+3]=(1-2*(param.q[2]*param.q[2]+param.q[3]*param.q[3]))*mag[0]+2*(param.q[1]*param.q[2]+param.q[0]*param.q[3])*mag[1]+2*(param.q[1]*param.q[3]-param.q[0]*param.q[2])*mag[2];
     param.H[j+12]=-2*mag[2]*param.q[2];
     break;

    case 1:
     param.h[j]=2*(param.q[2]*param.q[3]+param.q[0]*param.q[1]);
     param.H[j]=2*param.q[3];
     param.h[j+3]=(2*(param.q[1]*param.q[2]-param.q[0]*param.q[3])*mag[0]+(1-2*(param.q[1]*param.q[1]+param.q[3]*param.q[3]))*mag[1]+2*(param.q[2]*param.q[3]+param.q[0]*param.q[1])*mag[2]);
     param.H[j+12]=2*mag[2]*param.q[3]; 
     break;

    case 2:
     param.h[j]=1-2*(param.q[1]*param.q[1]+param.q[2]*param.q[2]);
     param.H[j]=-2*param.q[0];
     param.h[j+3]=(2*(param.q[1]*param.q[3]+param.q[0]*param.q[2])*mag[0]+2*(param.q[2]*param.q[3]-param.q[0]*param.q[1])*mag[1]+(1-2*(param.q[1]*param.q[1]+param.q[2]*param.q[2]))*mag[2]);
     param.H[j+12]=(-4)*mag[0]*param.q[2]-2*mag[2]*param.q[0];
     break; 
 
    case 3:
     param.H[j]=2*param.q[1];
     param.H[j+12]=(-4)*mag[0]*param.q[3]+2*mag[2]*param.q[1]; 
     break;

    case 4:
     param.H[j]=2*param.q[1];
     param.H[j+12]=2*(-1*mag[0]*param.q[3]+mag[2]*param.q[1]); 
     break;

    case 5:
     param.H[j]=2*param.q[0];
     param.H[j+12]=2*(mag[0]*param.q[2]+mag[2]*param.q[0]); 
     break;

    case 6:
     param.H[j]=2*param.q[3];
     param.H[j+12]=2*(mag[0]*param.q[1]+mag[2]*param.q[3]); 
     break;

    case 7:
     param.H[j]=2*param.q[2];
     param.H[j+12]=2*(-1*mag[0]*param.q[0]+mag[2]*param.q[2]); 
     break;

    case 8:
     param.H[j]=2*param.q[0];
     param.H[j+12]=2*mag[0]*param.q[2];
     break;

    case 9:
     param.H[j]=-2*param.q[1];
     param.H[j+12]=2*(mag[0]*param.q[3]-2*mag[2]*param.q[1]); 
     break;

    case 10:
     param.H[j]=-2*param.q[2];
     param.H[j+12]=2*(mag[0]*param.q[0]-2*mag[2]*param.q[2]); 
     break;

    case 11:
     param.H[j]=2*param.q[3];
     param.H[j+12]=2*mag[0]*param.q[1];
     break;
    }
  }
 } 
 void correction(void){
  // measurement residual (have to set h and H first)
   mat_sub(param.z,param.h,param.v,6);
  // Serial.println(" v");
  // print_mat(param.v,1,6);

  // measurement covariance
   mat_trans(param.H,param.temp,6,4);
   //Serial.println("H");
   //print_mat(param.H,6,4);
   //Serial.println("H^T");
   //print_mat(param.temp,4,6);
   
   mat_multiply(param.H,param.P,param.temp2,6,4,4);    //HP
   //Serial.println("HP");
   //print_mat(param.temp2,4,6);
   
   mat_multiply(param.temp2,param.temp,param.temp3,6,4,6);     //HPH^T     
   //Serial.println("HPH^T");
   //print_mat(param.temp3,6,6);
   mat_add(param.temp3,param.R,param.S,36);           //HPH^T+R
   //Serial.println("HPH^T+R");
   //print_mat(param.S,6,6);

  // kalman gain
   mat_inverse(param.S,param.S_inv,6);   //  need to be edit later

  //Serial.println(" S_inv");
  //print_mat(param.S_inv,6,6);
   
   mat_multiply(param.P, param.temp,param.temp2,4,4,6);    // PH^T
   mat_multiply(param.temp2,param.S_inv,param.K,4,6,6);    //K=PH^TS^-1 
  // state update
   mat_multiply(param.K,param.v,param.temp,4,6,1);        //kv
  //Serial.println(" Kv");
  //print_mat(param.K,1,4);
  mat_add(param.q,param.temp,param.q,4);      //q+kv

  // covariance update
   mat_multiply(param.K,param.H,param.temp,4,6,4);    //KH
   //Serial.println(" KH");
   //print_mat(param.temp,4,4);
 
   mat_sub(I_4by4,param.temp,param.temp,16);         //I-KH
   //Serial.println(" I-KH");
   //print_mat(param.temp,4,4);

   mat_multiply(param.temp,param.P,param.temp3,4,4,4);   //(I-KH)*P
   mat_copy(param.temp3,param.P,16);

   //Serial.println(" Next P");
   //print_mat(param.P,4,4);
   
  FTYPE buffer=InvSqrt(param.q[0]*param.q[0]+param.q[1]*param.q[1]+param.q[2]*param.q[2]+param.q[3]*param.q[3]);
  param.q[0]=param.q[0]*buffer;
  param.q[1]=param.q[1]*buffer;
  param.q[2]=param.q[2]*buffer;
  param.q[3]=param.q[3]*buffer;
 }   // combine measure covariance,kalman gain, state update, covariance update
    
};


