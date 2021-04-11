#ifndef SOURCE
#define SOURCE
#include "matrix.h"
#include <QString>
int  hoare(double r[], int l, int h);
void quicksort(double clock[], int epoch) ;  //快速排序算法？？？没细看

int  FixedClock(double clock[], int epoch) ; //钟差修复求斜率线性拟合

void  ClearGross( int FitNum, double phase[], int m_gross, int TempGrossNum); //对相位、频率采用进行MAD中位数法进行粗差探测
double Accuracyy(double t[], double phase[], int N); //频率准确度

double QuadraticPolynomial(double t[], double phase[], int FitNum);      //序贯平差进行二次多项式参数估计求解频率漂移率

int TotalHadamard(double phase[],int m_interval,int m_time,int FitNum,double &TempStable);     //哈达玛总方差函数声明

int ConvertWeeksec2Daysec(const int &Wsec);      //周内秒转换天内秒

void SetPath(int& year,int& month,int& day,int& hour,QString& filepath); //设置文件路径

int ReadData(const QString& filepath,double **data);            //读取文件 将钟差存入data


#endif // SOURCE

