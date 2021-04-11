#include "source.h"
#include <QFile>
#include <QString>
#include <QTextStream>

/*快速排序主体框架算法
input： double clock[]  钟差数据
        int    epoch    历元数目
output：clock[]
*/
void quicksort(double clock[], int epoch)
{
    int j, h, i;
    j = 1;
    h = epoch;   //历元数
    int tag = 1;
    int top = 0;
    int *s0;
    int *s1;
    s0 = new int[epoch];
    s1 = new int[epoch];
    do
    {
        while (j<h)
        {
            i = hoare(clock, j, h);
            top++;
            s0[top] = i + 1;
            s1[top] = h;
            h = i - 1;
        }
        if (top == 0)
        {
            tag = 0;
        }
        else
        {
            j = s0[top];
            h = s1[top];
            top--;
        }
    } while (tag == 1);
    delete[]s0;
    delete[]s1;
}
/**********快速排序分区处理函数**********
input:  double r[]
        int l
        int h
output：
*/
int  hoare(double r[], int l, int h)
{
    int i, j;
    double x;
    i = l;
    j = h;
    x = r[i - 1];
    do
    {
        while ((i < j) && (r[j - 1] >= x))
        {
            j--;
        }
        if (i < j)
        {
            r[i - 1] = r[j - 1];
            i++;
        }
        while ((i < j) && (r[i - 1] <= x))
        {
            i++;
        }
        if (i < j)
        {
            r[j - 1] = r[i - 1];
            j--;
        }
    } while (i<j);
    r[i - 1] = x;
    return(i);
}

/***********一次多项式修复无效钟差函数**********
input： double clock[] 某颗卫星钟差
        int epoch      历元数目
output：clock[]
*/
int  FixedClock(double clock[], int epoch)   //利用一次多项式进行钟差拟合，修复无效钟差数据
{
    int i = 0;
    int sum = 1;
    int  *usedate;                //使用的日期
    usedate = new int[epoch];     //历元数目
    for (int j = 0; j<epoch; j++) //初始化
    {
        usedate[j] = 0;
    }
    for (i = 1; i<epoch; i++)
    {
        if (fabs(clock[i])>1e-30)//如果是有效钟差，usedata[ ]=i等于数组号.若是无效钟差，usedate[ ]=0
        {
            usedate[sum] = i;
            sum++;
        }
    }
    for (i = 1; i<sum + 1; i++)
    {
        if (usedate[i] - usedate[i - 1] > 1)  //如果相邻历元相减大于1，认为该历元的钟差含有无效数据
        {
            int elem0 = usedate[i - 1];//前一个历元的数组号
            int elem1 = usedate[i];    //该历元的数组号
            double b00;
            b00 = (clock[elem1] - clock[elem0]) / (elem1 - elem0);//斜率
            for (int k = elem0 + 1; k<elem1; k++)
            {
                clock[k] = clock[elem0] + b00*(k - elem0);//
            }
        }
    }
    delete[]usedate;
    return 1;
}

/**********粗差探测*********
input:int FitNum         拟合历元数目
      double phase[]     钟差
      int m_gross        3
      int TempGrossNum   临时粗差数目数组
*/
void  ClearGross( int FitNum, double phase[], int m_gross, int TempGrossNum)
{
    //对相位数据进行粗差探测
    int i;
    //临时相位数组
    double *TempPhase;
    TempPhase = new double [FitNum];
    for (i = 0; i < FitNum; i++)
    {
        TempPhase[i] = phase[i];//相位值给临时相位数组
    }
    //相位序列的中间数
    double m;
    //快速排序算法
    quicksort(TempPhase, FitNum);
    m = TempPhase[FitNum/2];//钟差数据中位数
    //中位数
    double MAD=0.0;
    for (i = 0; i < FitNum; i++)
    {
        TempPhase[i] = fabs(phase[i] - m) / 0.6745;
    }
    quicksort(TempPhase, FitNum);
    MAD = TempPhase[FitNum/2];
    //探测粗差
    TempGrossNum = 0;
    for( i=0; i<FitNum; i++)
    {
        if(fabs(phase[i]-m) > m_gross*MAD)//粗差，m_gross=3；表示3倍MAD为限差
        {
            phase[i] = 0.0;
            TempGrossNum++;
        }
    }
    FixedClock(phase, FitNum);//修复
    //对频率数据进行粗差探测
    //临时频率数组
    double *TempFre;
    double *TempFre1;
    TempFre1 = new double [FitNum];
    TempFre = new double [FitNum];
    for( i=0; i<FitNum-1; i++)           //获取频率数据
    {
        TempFre1[i] = phase[i+1]-phase[i];
        TempFre[i] = phase[i+1]-phase[i];
    }
    //频率序列的中间数
    double fm=0.0;
    quicksort(TempFre1, FitNum-1);
    fm = TempFre1[(FitNum-1)/2];
    //中位数
    double *TempFreMad;
    TempFreMad = new double [FitNum];
    double FMAD=0.0;
    for (i = 0; i < FitNum - 1; i++)
    {
        TempFreMad[i] = fabs(TempFre[i] - fm) / 0.6745;
    }
    quicksort(TempFreMad, FitNum-1);
    FMAD = TempFreMad[(FitNum-1)/2];
    //探测粗差  严格的方法  大部分粗差能被探测
    for( i=0; i<FitNum-2; i++)
    {
        if((fabs(TempFre[i]-fm) > m_gross*FMAD)&&(fabs(TempFre[i+1]-fm) > m_gross*FMAD))//粗差
        {
            phase[i+1] = 0.0;
            TempGrossNum++;
        }
    }
    //开始两个数据的判断
    if((fabs(TempFre[0]-fm) > m_gross*FMAD)&&(fabs(TempFre[1]-fm) < m_gross*FMAD))//粗差
    {
        phase[0] = 0.0;
        TempGrossNum++;
    }
    //结尾两个数据的判断
    if((fabs(TempFre[FitNum-2]-fm) > m_gross*FMAD)&&(fabs(TempFre[FitNum-3]-fm) < m_gross*FMAD))//粗差
    {
        phase[FitNum-1] = 0.0;
        TempGrossNum++;
    }
    FixedClock(phase, FitNum);//修复
//释放动态数组
    delete []TempPhase;		TempPhase = NULL;
    delete []TempFre;		TempFre= NULL;
    delete []TempFre1;		TempFre1= NULL;
    delete []TempFreMad;	TempFreMad= NULL;
}


/**********计算频率准确度*********
input: double t[]         时间序列
       double phase[]     钟差序列
       int N              拟合数目
output:Kt                 频率准确度
notes:计算方法参考：高为广《北斗系统在轨卫星钟性能评估方法及结论》论文
*/
double Accuracyy(double t[], double phase[], int N)
{
    double tt[290] = { 0.0 };    //天内秒数组及初始化
    double AverageX = 0.0;       //钟差均值及初始化
    double Averaget = 0.0;       //时间均值及初始化
    double sumX = 0.0;           //钟差求和及初始化
    double sumt = 0.0;           //时间求和及初始化
    double Kt = 0.0;             //频率准确度及初始化
    double numerator = 0.0;      //计算频率准确度公式的分子及初始化
    double Denominator = 0.0;    //计算频率准确度公式的分母及初始化
    for (int k = 0; k < N; k++)
    {
        tt[k] = t[k] * 300;       //时间序列转为天内秒
        sumX = sumX + phase[k];   //逐个历元，相位求和
        sumt = sumt + tt[k];      //逐个历元，时间序列求和
    }
    AverageX = sumX / N;      //相位序列均值
    Averaget = sumt / N;      //时间序列均值
    for (int kk = 0; kk < N; kk++)
    {
        numerator = numerator + (phase[kk] - AverageX)*(tt[kk] - Averaget);    //计算频率准确度公式的分子
        Denominator = Denominator + (tt[kk] - Averaget)*(tt[kk] - Averaget);   //计算频率准确度公式的分母
    }
    Kt = numerator / Denominator;                                              //计算频率准确度
    return Kt;                                                                 //返回频率准确度
}

/*****************序贯平差进行二次多项式参数估计*****************
input： double t[]         时刻t的历元
        double phase[]     钟差
        int FitNum         历元数
        double x2          频率漂移率
        double residual[]  拟合残差数组
output：Fitprecision
*/
double QuadraticPolynomial(double t[], double phase[], int FitNum)
{
    double x2 = 0;
    MAT X0;       //定义X0矩阵
    MAT Q0;       //定义Q0矩阵
    while (1)     //前3个历元求解
    {
        //V=Ax-L
        MAT A(3, 3);//定义A矩阵。3行3列
        A.Set0();   //A矩阵初始化
        A.Set1();   //A矩阵设置为1
        A.SetElem(0, 1, t[0]); A.SetElem(0, 2, t[0] * t[0] / 2.0);//A矩阵赋值
        A.SetElem(1, 1, t[1]); A.SetElem(1, 2, t[1] * t[1] / 2.0);
        A.SetElem(2, 1, t[2]); A.SetElem(2, 2, t[2] * t[2] / 2.0);
        MAT L(3, 1);               //定义L矩阵，3行1列
        L.SetElem(0, 0, phase[0]); //L矩阵赋值
        L.SetElem(1, 0, phase[1]);
        L.SetElem(2, 0, phase[2]);
        /*矩阵：
        V  =         B       x   -   L
        x1     1 t1 t1*t1/2  a0    钟差1
        x2  =  1 t2 t2*t2/2  a1  + 钟差2
        x3     1 t3 t3*t3/2  a2    钟差3
        */
        Q0 = inverse1new(T(A)*A); //初始协方差，权为1
        X0 = Q0*T(A)*L;           //前三个钟差数据的a0，a1，a2
        break;
    }
    for (int i = 3; i<FitNum; i++)//序贯平差，逐个历元平差
    {
        MAT Aa(1, 3);    //定义Aa矩阵，1行3列
        Aa.Set1();
        MAT Ll(1, 1);    //定义L1矩阵，1行1列
        Ll.Set0();
        MAT J;           //定义J矩阵
        Aa.SetElem(0, 0, 1.0);                  //Aa矩阵赋值
        Aa.SetElem(0, 1, t[i]);
        Aa.SetElem(0, 2, t[i] * t[i] / 2.0);
        Ll.SetElem(0, 0, phase[i]);             //L1矩阵赋值

        MAT NN = T(Aa)*Aa;               //NN=AT*A
        MAT PP = inverse1(Q0);           //Q0求逆
        MAT Qk = inverse1(PP + NN);      //Qk
        MAT Xk = Qk*(T(Aa)*Ll + PP*X0);  //Xk
        X0 = Xk;
        Q0 = Qk;
    }
    x2 = X0.GetElem(2,0);                       //频率漂移率赋给x2
    return x2;
}

/****************哈达玛总方差计算稳定度****************
input： double phase[]       钟差
        int    m_interval    采样间隔：300秒
        int    m_time        平滑因子
        int    FitNum        拟合数据
        double &TempStable   稳定度
output：double &TempStable   稳定度
*/
int TotalHadamard(double phase[],int m_interval,int m_time,int FitNum,double &TempStable)
{
    //用相位数据计算哈达玛总方差
    double interval = m_interval;//采样间隔
    int m = m_time;              //平滑因子 对话框输入
    int N=FitNum;                //拟合数据个数
    int i;
    int j;
    int k;
    double c2=0.0;
    double sta =0.0;
    double sta1=0.0;
    double *phase1;                   //phase1存放3m点子序列
    phase1=new double[3*m+1];
    double *phase2;
    phase2=new double[9*m+1];         //phade2存放9m点延伸序列
    k=int(3*m/2);
    if (N - 4 * m < 0)
    {
        return 0;
    }
    for(j=0;j<N-3*m;j++)
    {
        if (m != 1)
        {
            c2 = (phase[j] - phase[j + k] - phase[3 * m - k + j] + phase[3 * m + j]) / k / (3 * m - k);//趋势项C2（见公式）
        }
        else                                   //采用反推的方法或取得c2值
        {
            c2 = (phase[j + 3] - 2 * phase[j + 2] + phase[j + 1]) +
                (phase[j + 3] - 3 * phase[j + 2] + 3 * phase[j + 1] - phase[j])*(sqrt(3) - 1) / 2;//不理解？？？
        }

        for (i = 0; i < 3 * m + 1; i++)                //3m+1点子序列生成
        {
            phase1[i] = phase[i + j] - c2*i*(i - 3 * m) / 2;
        }
        for (i = 0; i < 3 * m; i++)                  //9m+1点子序列生成
        {
            phase2[3 * m - i - 1] = 2 * phase1[0] - phase1[i + 1];
        }
        for (i = 3 * m; i < 6 * m + 1; i++)
        {
            phase2[i] = phase1[i - 3 * m];
        }
        for (i = 6 * m + 1; i < 9 * m + 1; i++)
        {
            phase2[i] = 2 * phase1[3 * m] - phase1[3 * m - (i - 6 * m)];
        }
        for( i=0;i<6*m;i++)
        {
            sta+=(phase2[i+3*m]-3*phase2[i+2*m]+3*phase2[i+m]-phase2[i])*
                    (phase2[i+3*m]-3*phase2[i+2*m]+3*phase2[i+m]-phase2[i]);
        }
        sta=sta/6/m;
        sta1+=sta;
        sta=0;
    }
    sta1=sta1/6/(N-3*m)/m/interval/m/interval;
    TempStable = sqrt(sta1);
    sta1=0.0;
    delete []phase1; phase1=NULL;
    delete []phase2; phase2=NULL;
     return 1;
}
int ConvertWeeksec2Daysec(const int &Wsec)
{
    int Dsec = Wsec%86400;
    return Dsec;
}

void SetPath(int& year,int& month,int& day,int& hour,QString& filepath)
{
    QString Syear,Smonth,Sday,Shour,Smin;
    Syear = QString ::number(year,10);
    Smonth = QString("%1").arg(month,2,10,QChar('0'));
    Sday = QString("%1").arg(day,2,10,QChar('0'));
    Shour = QString("%1").arg(hour,2,10,QChar('0'));

    filepath =  Syear+Smonth+Sday+Shour;

}

int ReadData(const QString& filepath,double **data)
{
    QFile file(filepath);
    if(!file.open( QIODevice::ReadOnly))
    {
       return 0;
    }
    int flag[36]={0};
    QTextStream read_clock(&file);
    while (!read_clock.atEnd())
    {
        QString str = "";
        str = read_clock.readLine();
        if(str.mid(0,1)=="4")
        {
            int count = str.mid(10,10).toInt();
            int prn = str.mid(9,1).toInt();
            if(flag[prn]==0)
            {
                flag[prn] = count;
                for(int i=0;i<count;i++)
                {
                    QString str = ""; int D_sec;
                    str = read_clock.readLine();
                   // int BD_week = str.mid(0,3).toInt();
                    int BD_sec = str.mid(4,8).toInt();
                    double clock = str.mid(13,20).toInt();
                    D_sec = ConvertWeeksec2Daysec(BD_sec);
                    data[prn][D_sec-1] = clock;
                }
            }
            else if(count>flag[prn])
            {
                flag[prn] = count;
                for(int i=0;i<count;i++)
                {
                    QString str = ""; int D_sec;
                    str = read_clock.readLine();
                   // int BD_week = str.mid(0,3).toInt();
                    int BD_sec = str.mid(4,8).toInt();
                    double clock = str.mid(13,20).toInt();
                    D_sec = ConvertWeeksec2Daysec(BD_sec);
                    data[prn][D_sec-1] = clock;
                }
            }else
            {
                for(int i=0;i<count;i++)
                {
                    read_clock.readLine();     //钟差已有跳过
                }
            }
        }

    }
    file.close();
    return 1;
}

