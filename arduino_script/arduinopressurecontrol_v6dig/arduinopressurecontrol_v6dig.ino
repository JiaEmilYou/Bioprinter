/*control the pressure controller and conmunicate with the PC
 * Get the pressure measurement from an external pressure
 */
const int analogInPin=A0;
const int analogOutPin=3; 
const int ExternalPressureIn=A3;
unsigned long time1;
int Vapply;
int Vmeasure;
int VEmeasure = 2;//(External measurement)
String printout;
int psensor=120;
//byte psensor=0xF1;
//byte psensor=0x78;
#include <Wire.h>


void setup() 
{
  // put your setup code here, to run once:
  pinMode(analogInPin,INPUT);
  pinMode(analogOutPin,OUTPUT);
  Serial.begin(115200);
  Wire.begin(); // join i2c bus (address optional for master)
  //digitalWrite(A4,HIGH);//pullup resistor
  //digitalWrite(A5,HIGH);//pullup resistor
}

void loop() 
{
  if (Serial.available() > 0) 
  {
    Vapply=Serial.read();
    analogWrite(analogOutPin,Vapply); 
  }
  Vmeasure=analogRead(analogInPin);
  Vmeasure=Vmeasure*32;
  Wire.beginTransmission(psensor);
  Wire.requestFrom(psensor,2);
  while (Wire.available()) 
  {
    //char c = Wire.read();
    //Serial.print(c);
    //Serial.println("Read from Pressure Sensor!");
    VEmeasure=(Wire.read()<<8)|(Wire.read());
    //Serial.println(VEmeasure);
  }
  Wire.endTransmission();
  /*
  //delay();
   */
  //VEmeasure=analogRead(ExternalPressureIn);
 
  
  time1 = millis();//updating time
  printout = "t="+String(time1)+" p="+VEmeasure;
  //printout=int(VEmeasure);
  Serial.println(printout);
  delay(5);
  
}
/*
unsigned long time;
void setup(){
  Serial.begin(9600);
}

void loop(){
  Serial.print("time=");
  time = millis();
  Serial.println(time);    //prints time since program started
  delay(500);
}

*/
