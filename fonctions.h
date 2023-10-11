#ifndef fonctions_h
#define fonctions_h


/************************
 fonctions conditions 
 aux limites et exactes
*************************/

float T_ex(float x, float y);
float f_0(float x, int m, float alpha);
float f_3_ex(float x);
float q_0(float x);

/************************
 fonctions a approcher
*************************/

float inte_h(float x, int m, float alpha);
float h(float x);
float inte_f3(float x, int m, float alpha);
float f_3(float x, float alpha);

/*************************
 fonctions a integrer pour
 les coefficients
**************************/

float B0_f3_f0(float x, int m, float alpha);
float Am_f3_f0(float x, int m, float alpha);
float Bm_f3_f0(float x, int m, float alpha);

/************************
 fonctions coefficients
*************************/

float A_0();
float B_0(float alpha);
float A_m(int m, float alpha);
float B_m(int m, float alpha);

/************************
 fonction T tilde approchee
*************************/

float T_tilde(float x, float y, float alpha);



#endif