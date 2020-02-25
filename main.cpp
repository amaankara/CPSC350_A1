#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <ctype.h>

using namespace std;
//Declaring global variables
int LineNum = 0;                        //calculating number of lines
int charNum = 0;                        //calculating number of characters
int dnaCount = 0;                       //calculating number of DNA strings
int Acnt = 0;                           //calculating number of As
int Ccnt = 0;                           //calculating number of Cs
int Tcnt = 0;                           //calculating number of Ts
int Gcnt = 0;                           //calculating number of Gs

float mean =0.0;                           //calculating number of mean
float var = 0.0;                           //calculating number of variance
float sd =0.0;                           //calculating number of standard deviation

//Mean of A,C,T,G
float meanA = 0.0;
float meanC = 0.0;
float meanT= 0.0;
float meanG = 0.0;


//count of all possible DNA string
int AA = 0;
int AC = 0;
int AT = 0;
int AG = 0;
int CA = 0;
int CC = 0;
int CT = 0;
int CG = 0;
int TA = 0;
int TC = 0;
int TT = 0;
int TG = 0;
int GA = 0;
int GC = 0;
int GT = 0;
int GG = 0;

//mean of each DNA string
float meanAA = 0.0;
float meanAC = 0.0;
float meanAT= 0.0;
float meanAG = 0.0;
float meanCA = 0.0;
float meanCC = 0.0;
float meanCT= 0.0;
float meanCG = 0.0;
float meanGA = 0.0;
float meanGC = 0.0;
float meanGT= 0.0;
float meanGG = 0.0;
float meanTA = 0.0;
float meanTC = 0.0;
float meanTT= 0.0;
float meanTG = 0.0;

//counting the DNA
void findString(string s){
  charNum += s.length();
  int num =0;
  num= s.length();
  //counting the amount of A,C,T and Gs
  for(int i =0; i< num;++i){
    char c = s[i];
    //cout << charNum << endl;
    if(c =='A'){
      Acnt++;
    }
    else if(c=='C'){
      Ccnt++;
    }
    else if(c=='T'){
      Tcnt++;
    }
    else if(c=='G'){
      Gcnt++;
    }
  }
  //counting all possible DNA stings
  for (int i = 1; i < num; ++i)
  {
    char first = s[i-1];
    char second = s[i];
    if((first =='A') && (second=='A')){
      AA++;
      dnaCount++;
    }
    else if((first =='A') && (second=='C')){
      AC++;
      dnaCount++;
    }
    else if((first =='A') && (second=='G')){
      AG++;
      dnaCount++;
    }
    else if((first =='A') && (second=='T')){
      AT++;
      dnaCount++;
    }
    else if((first =='C') && (second=='A')){
      CA++;
      dnaCount++;
    }
    else if((first =='C') && (second=='C')){
      CC++;
      dnaCount++;
    }
    else if((first =='C') && (second=='G')){
      CG++;
      dnaCount++;
    }
    else if((first ='C') && (second='T')){
      CT++;
      dnaCount++;
    }
    else if((first ='G') && (second='A')){
      GA++;
      dnaCount++;
    }
    else if((first ='G') && (second='C')){
      GC++;
      dnaCount++;
    }
    else if((first ='G') && (second='G')){
      GG++;
      dnaCount++;
    }
    else if((first ='G') && (second='T')){
      GT++;
      dnaCount++;
    }
    else if((first ='T') && (second='A')){
      TA++;
      dnaCount++;
    }
    else if((first ='T') && (second='C')){
      TC++;
      dnaCount++;
    }
    else if((first ='T') && (second='G')){
      TG++;
      dnaCount++;
    }
    else if((first ='T') && (second='T')){
      TT++;
      dnaCount++;
    }
  }
  //adding the number of lines
  LineNum ++;
  //cout<<LineNum<<endl;
}

//calculating the probability of each DNA
void probability(){
  mean = (float)charNum/(float)LineNum;
  //cout<<mean<<endl;
  var = (((float)charNum*charNum) - ((float)mean*mean))/LineNum;
  sd = sqrt(var);

  //calcuting probability of A,C,T and G
  meanA = (float)Acnt/(float)charNum;
  meanC = (float)Ccnt/(float)charNum;
  meanG = (float)Gcnt/(float)charNum;
  meanT = (float)Tcnt/(float)charNum;

  //calcuting the probability of each possible DNA string
  meanAA = (float)AA/(float)dnaCount;
  meanAC = (float)AC/(float)dnaCount;
  meanAT = (float)AT/(float)dnaCount;
  meanAG = (float)AG/(float)dnaCount;
  meanCA = (float)CA/(float)dnaCount;
  meanCC = (float)CC/(float)dnaCount;
  meanCT = (float)CT/(float)dnaCount;
  meanCG = (float)CG/(float)dnaCount;
  meanTA = (float)TA/(float)dnaCount;
  meanTC = (float)TC/(float)dnaCount;
  meanTT = (float)TT/(float)dnaCount;
  meanTG = (float)TG/(float)dnaCount;
  meanGA = (float)GA/(float)dnaCount;
  meanGC = (float)GC/(float)dnaCount;
  meanGT = (float)GT/(float)dnaCount;
  meanGG = (float)GG/(float)dnaCount;
}

//calcuting the Gaussian variable
float gauDistribution(float avg, float standev){
  //random variable a and b between [0,1)
  float A_rand = (float)rand()/(float)RAND_MAX;
  float B_rand = (float)rand()/(float)RAND_MAX;
  //cout<< A_rand << " " << B_rand <<endl;
  float lg = log(A_rand);
  float C = sqrt((-2)*lg*cos((M_PI)*(2)*B_rand));
  float D = standev*C + avg;
  //if d < 0 the making it positive
  if(D<=0){
    D = (-1)*D;
  }

  return D;
}

void writeToFile(){
  //calling the probability() function
  probability();
  //outputing results to the file
  ofstream outFile("AmaanKarachiwala.out");
  outFile << "Amaan Karachiwala; 2326810; A1" << endl;
  outFile << "The mean number of lines is " << mean << endl;
  outFile << "The standard deviation is " << sd << endl;
  outFile << "The probability of A is " << Acnt << endl;
  outFile << "The probability of C is  " << meanC << endl;
  outFile << "The probability of T is " << meanT << endl;
  outFile << "The probability of G is  " << meanG << endl;

  outFile << "The probability of AA is " << meanAA << endl;
  outFile << "The probability of AC is " << meanAC << endl;
  outFile << "The probability of AT is " << meanAT << endl;
  outFile << "The probability of AG is " << meanAG << endl;
  outFile << "The probability of CA is  " << meanCA << endl;
  outFile << "The probability of CC is  " << meanCC << endl;
  outFile << "The probability of CT is  " << meanCT << endl;
  outFile << "The probability of CG is  " << meanCG << endl;
  outFile << "The probability of TA is " << meanTA << endl;
  outFile << "The probability of TC is " << meanTC << endl;
  outFile << "The probability of TT is " << meanTT << endl;
  outFile << "The probability of TG is " << meanTG << endl;
  outFile << "The probability of GA is  " << meanGA << endl;
  outFile << "The probability of GC is  " << meanGC << endl;
  outFile << "The probability of GT is  " << meanGT << endl;
  outFile << "The probability of GG is  " << meanGG << endl;
  string writeline = "";
  int loopCount = 0;
  int a = 0;
  int c = 0;
  int t = 0;
  int g = 0;

  //creating the new strings
  for(int i=0; i<1000; ++i){
    //creating the Gaussian varibale
    float gaussainVar = gauDistribution(mean, sd);
    int writeA = (int)((float)Acnt/charNum)*gaussainVar;
    //cout<<sd  << " " << mean << " " << meanA << " " << gauDistribution << endl;
    int writeC = (int)((float)Ccnt/charNum)*gaussainVar;
    int writeT = (int)((float)Tcnt/charNum)*gaussainVar;
    int writeG = (int)((float)Gcnt/charNum)*gaussainVar;

    // if(writeA >= writeC && writeA >= writeT && writeA >= writeG){
    //   loopCount = writeA;
    // }
    // else if(writeC >= writeA && writeC >= writeT && writeC >= writeG){
    //   loopCount = writeC;
    // }
    // else if(writeT >= writeC && writeT >= writeA && writeT >= writeG){
    //   loopCount = writeT;
    // }
    //
    // else {
    //   loopCount = writeG;
    // }

    //adding A,C,T and G to the string
    for (int j = 0; j < writeA; ++j) {
      writeline+="A";
    }
    for (int j = 0; j < writeC; ++j) {
      writeline+="C";
    }
    for (int j = 0; j < writeT; ++j) {
      writeline+="T";
    }
    for (int j = 0; j < writeG; ++j) {
      writeline+="G";
    }

    // for (int j = 0; j < loopCount; ++j) {
    //   if(a < writeA){
    //     writeline += "A";
    //     a++;
    //   }
    //  if(c < writeC){
    //    writeline += "C";
    //    c++;
    //  }
    //  if(t < writeT){
    //    writeline += "T";
    //    t++;
    //  }
    //  if(g < writeG){
    //    writeline += "G";
    //    g++;
    //  }
    // }
    //cout<<writeline<<endl;

    //outputing to the file
    outFile << writeline << endl;
    //clearing the string
    writeline = "";
  }
}


int main(int argc, char **argv){
  //varibale to check if user is done or not
  char yesNo = 'y';
  //declaring local varibales
  string ln;

  //getting filename from user
  string filename = "";
  //running a while loop until the user is done
  while(yesNo=='y' || yesNo == 'Y'){
    //inputing the file from the user
    cout << "Enter the filename " << endl;
    cin >> filename;

    //opening the file line by line and sending it to another function to count DNA
    ifstream myfile(filename);
    while(getline(myfile,ln)){
      //cout << ln << endl;
      findString(ln);
    }
    //displayProb();
    //calling the function writeToFile();
    writeToFile();
    //asking user if they want to continue
    cout<<"Would you like to continue (Y/N)" << endl;
    cin >> yesNo;
  }
  return 0;
}
