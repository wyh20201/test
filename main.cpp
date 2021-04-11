#include <QString>
#include <QTextStream>
#include <QFile>
#include "source.h"
#include "stdio.h"

using namespace std;

int main()
{

    FILE *file = fopen("*.txt","r");
    QString Syear,Smonth,Sday,Shour,Smin;
    Syear = QString ::number(2019,10);
    Sday = QString("%1").arg(1,2,10,QChar('0'));
    Shour = Syear+Sday;
    return 0;
}

