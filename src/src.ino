#include <FaBo9Axis_MPU9250.h>
#include "ekf.h"

ekf orient;
FaBo9Axis fabo_9axis;

float ax,ay,az,gx,gy,gz,mx,my,mz;

FTYPE accel_measure[3];
FTYPE accelbias[3];

FTYPE mag_measure[3];
FTYPE magbias[3]={7.68,-9.406,13.154};   // the already calibration measurement value
FTYPE A_inv[9]={0.02765,0.00207,-0.0002464,0.00207,0.0285,-0.00064,-0.0002464,-0.00064,0.0285}; 

FTYPE gyro_measure[3];
FTYPE gyrobias[3];

FTYPE last_time=0;
FTYPE now,dt;

void setup() {
  Serial.begin(115200);
  Serial.println("RESET");
  Serial.println();

  Serial.println("configuring device.");
  fabo_9axis.configMPU9250(MPU9250_GFS_500,MPU9250_AFS_4G);
  if (fabo_9axis.begin()) {
    Serial.println("configured FaBo 9Axis I2C Brick");
  } else {
    Serial.println("device error");
    while(1);
  }
  
  // ekf setting the magnetometer is currectly usable( but not accurate enough), it still need improvement
  Serial.println("please place the sensor in flat");
  delay(1000);
  fabo_9axis.readMagnetXYZ(&mag_measure[1],&mag_measure[0],&mag_measure[2]);  // due to definition of direction for MPU9250 y->x x->y
  fabo_9axis.readAccelXYZ(&accel_measure[0],&accel_measure[1],&accel_measure[2]);
 // mat_sub(mag_measure,magbias,mag_measure,3);
 // mat_multiply(A_inv,mag_measure,mag_measure,3,3,1);
  acc_gyro_calib();
  mat_sub(accel_measure,accelbias,accel_measure,3);

  Serial.print("before calibration: ");
  Serial.print("mx: ");
  Serial.print(mag_measure[0]);
  Serial.print(" my: ");
  Serial.print(mag_measure[1]);
  Serial.print(" mz: ");
  Serial.println(mag_measure[2]);
 
  initial_attitude();

  // noise setting
  //gyro: 0.1 deg/s-rms=3.0462e-06 rad/s->square 
  // acc: 8mg/s^2 rms=0.0784 m/s^2-rms->square mag=0.6 muT/LSB->square

  orient.ekf_measumement_noise(0.0061,0.0061,0.0061,0.36,0.36,0.36); 
  orient.ekf_gyro_noise(0.000003,0.000003,0.000003); 

  
  

}

void loop(){
  fabo_9axis.readMagnetXYZ(&mag_measure[1],&mag_measure[0],&mag_measure[2]);  // mx=mz mz=mx 
  fabo_9axis.readAccelXYZ(&accel_measure[0],&accel_measure[1],&accel_measure[2]);
  fabo_9axis.readGyroXYZ(&gyro_measure[0],&gyro_measure[1],&gyro_measure[2]);
  //mat_sub(mag_measure,magbias,mag_measure,3);
  //mat_multiply(A_inv,mag_measure,mag_measure,3,3,1);

  now = (FTYPE) millis()/1000;
  dt = now - last_time;
  last_time = now;

  mat_sub(accel_measure,accelbias,accel_measure,3);
  mat_sub(gyro_measure,gyrobias,gyro_measure,3);

  // from deg/s scale to rad/s due to MPU9250 angular velocity direction, the wy wz has to be negative
  gyro_measure[0]=gyro_measure[0]*PI/180;
  gyro_measure[1]=-1*gyro_measure[1]*PI/180;   
  gyro_measure[2]=-1*gyro_measure[2]*PI/180;

  //Serial.print("prediction: ");
  orient.ekf_prediction(gyro_measure[0],gyro_measure[1],gyro_measure[2],dt);
  //Serial.print("correction: ");
  orient.ekf_correction(accel_measure,mag_measure);
  orient.ekf_qtn_to_Euler();

  Serial.print("Roll, Pitch, Yaw: ");
  Serial.print(orient.roll, 2);
  Serial.print(" ");
  Serial.print(orient.pitch, 2);
  Serial.print(" ");
  Serial.print(orient.yaw, 2);
  Serial.println("  "); 

}


void initial_attitude(void){
    orient.ekf_init(accel_measure[0],accel_measure[1],accel_measure[2],mag_measure[0],mag_measure[1],mag_measure[2]);
}

void acc_gyro_calib(void){
  //assume the sensor is horizontal
  Serial.println("accel gyro calib, please place it on horizontal plane");
  delay(200);
  fabo_9axis.readAccelXYZ(&accel_measure[0],&accel_measure[1],&accel_measure[2]);
  fabo_9axis.readGyroXYZ(&gyro_measure[0],&gyro_measure[1],&gyro_measure[2]);
  accelbias[0]=accel_measure[0];
  accelbias[1]=accel_measure[1];
  accelbias[2]=accel_measure[2]-1;
  gyrobias[0]=gyro_measure[0];
  gyrobias[1]=gyro_measure[1];
  gyrobias[2]=gyro_measure[2];
  Serial.println("accel gyro calib finish");
 
}

/*
  Serial.print("before calibration: ");
  Serial.print("mx: ");
  Serial.print(mag_measure[2]);
  Serial.print(" my: ");
  Serial.print(mag_measure[1]);
  Serial.print(" mz: ");
  Serial.println(mag_measure[0]);
 
  //  magnetometer calibration
  mat_sub(mag_measure,magbias,mag_measure,3);
  mat_multiply(A_inv,mag_measure,mag_measure,3,3,1);
  Serial.print("after calibration ");
  Serial.print("mx: ");
  Serial.print(mag_measure[2]);
  Serial.print(" my: ");
  Serial.print(mag_measure[1]);
  Serial.print(" mz: ");
  Serial.println(mag_measure[0]);

  fabo_9axis.readAccelXYZ(&accel_measure[0],&accel_measure[1],&accel_measure[2]);
  Serial.print("before calibration: ");
  Serial.print("accx: ");
  Serial.print(accel_measure[0]);
  Serial.print(" accy: ");
  Serial.print(accel_measure[1]);
  Serial.print(" accz: ");
  Serial.println(accel_measure[2]);
  
  mat_sub(accel_measure,accelbias,accel_measure,3);
  Serial.print("after calibration: ");
  Serial.print("accx: ");
  Serial.print(accel_measure[0]);
  Serial.print(" accy: ");
  Serial.print(accel_measure[1]);
  Serial.print(" accz: ");
  Serial.println(accel_measure[2]);
  delay(100);


   // loop part
  


*/