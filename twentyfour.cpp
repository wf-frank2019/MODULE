#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

int main()
{
	int card[4];
	char word[20];
	char c;
	printf("******************************\n");
	printf("	  欢迎来到24点游戏\n");
	do
	{
		getcard(card);		//随机生成卡牌 
		printf("你所得到的4个随机卡牌号为(用L代替10)；\n");
		printcard(card);			//输出卡牌号 
		printf("请按任意键查看结果（用L代替10）\n");
		getchar();		//获取所按的键 
		printresult(card,word); 		//输出所有运算结果 
		printf("\n继续游戏请按回车键，结束请按其他键\n");
        c=getchar();//将所按的键的命令给C 
	}
	while(c=='\n');			//判断所按的键来决定是否继续游戏 
	return 0;

}

int getcard(int *card)						//产生随机数 
{
	int i;
	srand((unsigned)time(NULL));			 //利用系统时间作为种子产生随机数
    for(i=0;i<4;i++)
    {
        card[i]=rand()%12+1;			//用rand函数生成随机数，范围为1-13 
    }
    return;
}

int printcard(int *card)		//输出随机得到的卡牌号 
{
	int i;
	for(i=0;i<4;i++)
	{
	/*	if(card[i]==10)
		{
			printf("%d  ",card[i]);
		}
		else
		{*/
			printf("%c  ",getcardcard(card[i]));
		//}
	}
	printf("\n");
	return;
}

int getcardcard(int a)				//将随机数字转化为卡牌号 
{
	if(a==1)
	{
		return 'A';
	}
	else if(a<10)
	{
		return a+'0';
	}
	else if(a==10)
	{
		return 'L';
	 } 
	else if(a==11)
	{
		return 'J';
	}
	else if(a==12)
	{
		return 'Q';
	}
	else if(a==13)
	{
		return 'K';
	}
}



int printresult(int *card,char *word)			//得到计算结果
{	
	char b[4]={'+','-','*','/'};
	char operator[3];
	int i,j,k,t=0;
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			for(k=0;k<4;k++)
			{
				operator[0]=b[i];
				operator[1]=b[j];
				operator[2]=b[k];
				if(getresult(card,operator,word))
				t++;
			}
		}
	}
	if(t!=0)
	{
		printf("共有%d种解\n",t);
	}
	else
	{
		printf("该情况无解\n");
	}
	return;
}

void printword(int flag,int *card,char *operator,char *word)		//输出计算后的结果 
{
    char a=getcardcard(card[0]);
    char b=getcardcard(card[1]);
    char c=getcardcard(card[2]);
    char d=getcardcard(card[3]);
     
    switch(flag)
    {
        //1.((A*B)*C)*D
        case 1:
            printf("((%c%c%c)%c%c)%c%c=24\n",a,operator[0],b,operator[1],c,operator[2],d);
            break;
        //2.(A*(B*C))*D
        case 2:
            printf("(%c%c(%c%c%c))%c%c=24\n",a,operator[0],b,operator[1],c,operator[2],d);
            break;
        //3.(A*B)*(C*D)
        case 3:
            printf("(%c%c%c)%c(%c%c%c)=24\n",a,operator[0],b,operator[1],c,operator[2],d);
            break;
        //4.A*(B*(C*D))
        case 4:
            printf("%c%c(%c%c(%c%c%c))=24\n",a,operator[0],b,operator[1],c,operator[2],d);
            break;
        //5.A*((B*C)*D) 
        case 5:
            printf("%c%c((%c%c%c)%c%c)=24\n",a,operator[0],b,operator[1],c,operator[2],d);
            break;
        default:
            break;
    }
}

double getvalue(double num1,double num2,char operator)		//求两个数的计算结果 
{
    double result;
     
    switch(operator)
    {
        case '+':
        result=num1+num2;
        break;
        case '-':
        result=fabs(num1-num2);
        break;
        case '*':
        result=num1*num2;
        break;
        case '/':
        result=num1/num2;
        break;
        default :
        break;
    }
    return result;
}    
    
 
int getresult(int *card,char *operator,char *word)
{
    double t;
    //将计算值取到 
    int a=card[0]>10?1:card[0];	
    int b=card[1]>10?1:card[1];
    int c=card[2]>10?1:card[2];
    int d=card[3]>10?1:card[3];
     
    //穷举运算次序
    //1.((A*B)*C)*D
    t=0;
    t=getvalue(a,b,operator[0]);				//通过分层调用getvalue函数，达到括号的效果 
    t=getvalue(t,c,operator[1]);
    t=getvalue(t,d,operator[2]);
    
    if(fabs(t-24)<0.0001)		//fabs函数求绝对值 
    {
        printword(1,card,operator,word);
        return 1;
    }
    
    //2.(A*(B*C))*D
    t=0;
    t=getvalue(b,c,operator[1]);	//求得最内括号内的表达式计算结果 
    t=getvalue(a,t,operator[0]);	//再依次求得括号内的表达式计算结果 
    t=getvalue(t,d,operator[2]);	//求最外层表达式的结果 
    
    if(fabs(t-24)<0.0001)		//判断结果是否为24 
    {
        printword(2,card,operator,word);
        return 1;
    }
    
    //3.(A*B)*(C*D)
    t=0;
    t=getvalue(getvalue(a,b,operator[0]),getvalue(c,d,operator[2]),operator[1]);
    
    if(fabs(t-24)<0.0001)
    {
        printword(3,card,operator,word);
        return 1;
    }
    
    //4.A*(B*(C*D))
    t=0;
    t=getvalue(c,d,operator[2]);
    t=getvalue(b,t,operator[1]);
    t=getvalue(a,t,operator[0]);
    
    if(fabs(t-24)<0.0001)
    {
        printword(4,card,operator,word);
        return 1;
    }
    
    //5.A*((B*C)*D) 
    t=0;
    t=getvalue(b,c,operator[1]);
    t=getvalue(t,d,operator[2]);
    t=getvalue(a,t,operator[0]);
    
    if(fabs(t-24)<0.0001)
    {
        printword(5,card,operator,word);
        return 1;
    }
     return 0;
 }
