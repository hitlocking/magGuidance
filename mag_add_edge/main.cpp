#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/robust_kernel.h"
#include "g2o/core/robust_kernel_impl.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"

#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"

//#include <g2o/types/slam3d/vertex_se3.h>
//#include <g2o/types/slam3d/edge_se3.h>
#include <g2o/types/slam2d/types_slam2d.h>
// 使用 宏函数 声明边和顶点类型，注意注释掉了上面两个头文件 
//G2O_USE_TYPE_GROUP(slam3d);
//G2O_USE_TYPE_GROUP(slam2d); //2d平面

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

using namespace std;

using namespace g2o;

#define MAX_LINE 200
#define MAXITERATION 60	

void MyWrite(int num1, int num2, FILE *fpp);
float min3(float num1, float num2, float num3);
int min2(int num1, int num2);
void FileAnalysis(void);
void MissCheck(void);
//利用结构体存储数据
struct VERTEX	
{
	unsigned int id;
	float x;
	float y;
	float theta;
};
struct EDGE
{
	unsigned int fromid;
	unsigned int toid;
	float delta_x;
	float delta_y;
	float delta_t;
	float info_xx;	//信息矩阵
	float info_xy;
	float info_xt;
	float info_yy;
	float info_yt;
	float info_tt;
};
struct VERTEX vertex_all[4000];	//存储vertex数据的数组，4维，数组下标代表vertex的序号
struct EDGE edge_con[4000];		//存储edge数据的数组，11维，相邻边
struct EDGE edge_Ncon[30];		//存储edge数据的数组，11维，非相邻边
struct EDGE edge_wrong[50];		//存储edge数据的数组，11维，错误的edge测量
int NumOfedge_Ncon = 0;		//存储非相邻边的数量
int NumOfedge_Con = 0;		//存储相邻边的数量
int NumOfVertex = 0;		//存储节点的数量
int Cir1_Point1, Cir2_Point1, Cir3_Point1;
int Cir1_Point2, Cir2_Point2, Cir3_Point2;
int Cir1_Point3, Cir2_Point3, Cir3_Point3;
int Cir1_Point4, Cir2_Point4, Cir3_Point4;
int Cir1_Point5, Cir2_Point5, Cir3_Point5;


int main()
{
	FILE *fpp1;
	FILE *fpp2;
	char file_NS[] = "../data/input.g2o";	//fpp1
	char fi[] = "../data/added.g2o";					//fpp2
	char ch[256];
	if((fpp1=fopen(file_NS,"r")) == NULL){	//打开只读文件，该文件必须存在
		perror("fail to open input.g2o!");
		exit(1);
	}
	else printf("open succeed input.g2o!\n");
	
	if((fpp2=fopen(fi,"w+")) == NULL){	//打开可读写文件,删除文件全部内容，若不存在则建立文件
		perror("fail to open added.g2o!");
		exit(1);
	}
	else printf("open succeed added.g2o!\n");
	//把graph_3laps_NS.g2o中的数据复制到added.g2o中
	while(fgets(ch,256,fpp1)!=NULL){	//从stream中读取一行
		puts(ch);	//打印输出
		fputs(ch, fpp2);	//将字符串ch输出到文件fpp2中
	}		
	fclose(fpp1);
	FileAnalysis();
//以上完成了对所有数据的分析和提取，存储在三个结构体数组中，分别是：vertex_all，edge_Con，edge_Ncon
//下面提取了4个锁定点，也就是最初的约束边连接的点
	int LockPt[3][4];
	int pCnt;
	Cir1_Point1 = edge_Ncon[0].fromid;
	Cir2_Point1 = edge_Ncon[0].toid;
	for(pCnt=1; pCnt <= NumOfedge_Ncon; pCnt++){
		if(edge_Ncon[pCnt].fromid == edge_Ncon[0].toid){
			Cir3_Point1 = edge_Ncon[pCnt].toid;
		}
	}
	Cir1_Point2 = edge_Ncon[1].fromid;
	Cir2_Point2 = edge_Ncon[1].toid;
	for(pCnt=2; pCnt <= NumOfedge_Ncon; pCnt++){
		if(edge_Ncon[pCnt].fromid == edge_Ncon[1].toid){
			Cir3_Point2 = edge_Ncon[pCnt].toid;
		}
	}
	Cir1_Point3 = edge_Ncon[2].fromid;
	Cir2_Point3 = edge_Ncon[2].toid;
	for(pCnt=3; pCnt <= NumOfedge_Ncon; pCnt++){
		if(edge_Ncon[pCnt].fromid == edge_Ncon[2].toid){
			Cir3_Point3 = edge_Ncon[pCnt].toid;
		}
	}
	Cir1_Point4 = edge_Ncon[3].fromid;
	Cir2_Point4 = edge_Ncon[3].toid;
	for(pCnt=4; pCnt <= NumOfedge_Ncon; pCnt++){
		if(edge_Ncon[pCnt].fromid == edge_Ncon[3].toid){
			Cir3_Point4 = edge_Ncon[pCnt].toid;
		}
	}
	
	printf("point1:%d   %d   %d\n", Cir1_Point1, Cir2_Point1, Cir3_Point1);	
	printf("point2:%d   %d   %d\n", Cir1_Point2, Cir2_Point2, Cir3_Point2);
	printf("point3:%d   %d   %d\n", Cir1_Point3, Cir2_Point3, Cir3_Point3);
	printf("point4:%d   %d   %d\n", Cir1_Point4, Cir2_Point4, Cir3_Point4);
	LockPt[0][0] = Cir1_Point1; LockPt[1][0] = Cir2_Point1; LockPt[2][0] = Cir3_Point1;
	LockPt[0][1] = Cir1_Point2; LockPt[1][1] = Cir2_Point2; LockPt[2][1] = Cir3_Point2;
	LockPt[0][2] = Cir1_Point3; LockPt[1][2] = Cir2_Point3; LockPt[2][2] = Cir3_Point3;
	LockPt[0][3] = Cir1_Point4; LockPt[1][3] = Cir2_Point4; LockPt[2][3] = Cir3_Point4;
	int qCnt=0;

	int a;
	//对三圈中的第一、二、三、四段进行匹配
	for(a=0;a<3;a++){
	int i=0,j=0,k=0;
		while(i<(LockPt[0][a+1]-LockPt[0][a]) && j<(LockPt[1][a+1]-LockPt[1][a]) && k<(LockPt[2][a+1]-LockPt[2][a])){ 
			if( (fabs(edge_con[LockPt[0][a]+i].delta_x - edge_con[LockPt[1][a]+j].delta_x)<0.5) &&
				(fabs(edge_con[LockPt[1][a]+j].delta_x - edge_con[LockPt[2][a]+k].delta_x)<0.5) &&
				(fabs(edge_con[LockPt[0][a]+i].delta_x - edge_con[LockPt[2][a]+k].delta_x)<0.5) ){
				//如果第一圈的点和第二圈的点，以及第二圈的点和第三圈的点的差，都比较小，则这三个点是同一个点，用边连接起来。
				MyWrite(LockPt[0][a]+i, LockPt[1][a]+j, fpp2);
				MyWrite(LockPt[1][a]+j, LockPt[2][a]+k, fpp2);
				i++; j++; k++;
			}
			else{		//不插值，漏检的磁钉，在其他圈中把对应磁钉删除
				float temp;
				temp = min3( fabs(edge_con[LockPt[0][a]+i].delta_x), fabs(edge_con[LockPt[1][a]+j].delta_x), 
							 fabs(edge_con[LockPt[2][a]+k].delta_x) );
				if(temp == fabs(edge_con[LockPt[0][a]+i].delta_x)){
					edge_con[LockPt[0][a]+i+1].delta_x = edge_con[LockPt[0][a]+i+1].delta_x + edge_con[LockPt[0][a]+i].delta_x;
					printf("漏检磁钉：%d\n", LockPt[0][a]+i);
					i = i+1;
				}
				if(temp == fabs(edge_con[LockPt[1][a]+j].delta_x)){
					edge_con[LockPt[1][a]+j+1].delta_x = edge_con[LockPt[1][a]+j+1].delta_x + edge_con[LockPt[1][a]+j].delta_x;
					printf("漏检磁钉：%d\n", LockPt[1][a]+j);
					j = j+1;
				}
				if(temp == fabs(edge_con[LockPt[2][a]+k].delta_x)){
					edge_con[LockPt[2][a]+k+1].delta_x = edge_con[LockPt[2][a]+k+1].delta_x + edge_con[LockPt[2][a]+k].delta_x;
					printf("漏检磁钉：%d\n", LockPt[2][a]+k);
					k = k+1;
				}
			}
		}
	}
	//对0点到第一个锁定点之间的点进行匹配
	int i=LockPt[0][0],j=LockPt[1][0], k=LockPt[2][0];
	int NonUse1=0, NonUse2=0;
	while(i>=0){
		if( (fabs(edge_con[i].delta_x - edge_con[j].delta_x)<0.5) &&
			(fabs(edge_con[j].delta_x - edge_con[k].delta_x)<0.5) &&
			(fabs(edge_con[i].delta_x - edge_con[k].delta_x)<0.5) ){
				//如果第一圈的点和第二圈的点，以及第二圈的点和第三圈的点的差，都比较小，则这三个点是同一个点，用边连接起来。
				MyWrite(i, j, fpp2);
				MyWrite(j, k, fpp2);
				i--; j--; k--;
				NonUse1 = j; NonUse2 = k;
		}
		else{		//不插值，漏检的磁钉，在其他圈中把对应磁钉删除
			float temp;
			temp = min3( fabs(edge_con[i].delta_x), fabs(edge_con[j].delta_x), 
						 fabs(edge_con[k].delta_x) );
			if(temp == fabs(edge_con[i].delta_x)){
				edge_con[i-1].delta_x = edge_con[i-1].delta_x + edge_con[i].delta_x;
				printf("漏检磁钉：%d\n", i);
				i = i-1;
			}
			if(temp == fabs(edge_con[j].delta_x)){
				edge_con[j-1].delta_x = edge_con[j-1].delta_x + edge_con[j].delta_x;
				printf("漏检磁钉：%d\n", j);
				j = j-1;
			}
			if(temp == fabs(edge_con[k].delta_x)){
				edge_con[k-1].delta_x = edge_con[k-1].delta_x + edge_con[k].delta_x;
				printf("漏检磁钉：%d\n", k);
				k = k-1;
			}
		}
	}
	//对第五个锁定点到终点的点进行匹配
	i=LockPt[0][3],j=LockPt[1][3], k=LockPt[2][3];
	while(k<=NumOfedge_Con){
		if(k == NumOfedge_Con){
			MyWrite(j, k, fpp2);
			break;
		}	
		if( (fabs(edge_con[i].delta_x - edge_con[j].delta_x)<0.5) &&
			(fabs(edge_con[j].delta_x - edge_con[k].delta_x)<0.5) &&
			(fabs(edge_con[i].delta_x - edge_con[k].delta_x)<0.5) ){
				//如果第一圈的点和第二圈的点，以及第二圈的点和第三圈的点的差，都比较小，则这三个点是同一个点，用边连接起来。
				MyWrite(i, j, fpp2);
				MyWrite(j, k, fpp2);
				i++; j++; k++;
		}
		else{		//不插值，漏检的磁钉，在其他圈中把对应磁钉删除
			float temp;
			temp = min3( fabs(edge_con[i].delta_x), fabs(edge_con[j].delta_x), 
						 fabs(edge_con[k].delta_x) );
			if(temp == fabs(edge_con[i].delta_x)){
				edge_con[i+1].delta_x = edge_con[i+1].delta_x + edge_con[i].delta_x;
				printf("漏检磁钉：%d\n", i);
				i = i+1;
			}
			if(temp == fabs(edge_con[j].delta_x)){
				edge_con[j+1].delta_x = edge_con[j+1].delta_x + edge_con[j].delta_x;
				printf("漏检磁钉：%d\n", j);
				j = j+1;
			}
			if(temp == fabs(edge_con[k].delta_x)){
				edge_con[k+1].delta_x = edge_con[k+1].delta_x + edge_con[k].delta_x;
				printf("漏检磁钉：%d\n", k);
				k = k+1;
			}
		}
	}
	
	while(i<=NonUse1 && j<=NonUse2){
		if(fabs(edge_con[i].delta_x - edge_con[j].delta_x)<0.5){
			//如果第一圈的点和第二圈的点，以及第二圈的点和第三圈的点的差，都比较小，则这三个点是同一个点，用边连接起来。
			MyWrite(i, j, fpp2);		
			i++; j++;
		}
		else{		//不插值，漏检的磁钉，在其他圈中把对应磁钉删除
			float temp;
			temp = min2( fabs(edge_con[i].delta_x), fabs(edge_con[j].delta_x) );
			if(temp == fabs(edge_con[i].delta_x)){
				edge_con[i+1].delta_x = edge_con[i+1].delta_x + edge_con[i].delta_x;
				printf("漏检磁钉：%d\n", i);
				i = i+1;
			}
			if(temp == fabs(edge_con[j].delta_x)){
				edge_con[j+1].delta_x = edge_con[j+1].delta_x + edge_con[j].delta_x;
				printf("漏检磁钉：%d\n", j);
				j = j+1;
			}
		}
	}
				
	printf("done!\n");
	fclose(fpp1);
	fclose(fpp2);
	
	cout<< "Hello g2o"<<endl;
	// create the optimizer,构建求解器
	SparseOptimizer optimizer;	
	// create the linear solver， 构建线性求解器，CSparse
	BlockSolverX::LinearSolverType* linearSolver = new LinearSolverCSparse<BlockSolverX::PoseMatrixType>();
	// create the block solver on the top of the linear solver， 构建block slover
	BlockSolverX* blockSolver = new BlockSolverX(linearSolver);
	// create the algorithm to carry out the optimization， 构建算法，LM下降法
	OptimizationAlgorithmLevenberg* Algorithm = new OptimizationAlgorithmLevenberg(blockSolver);  

	// 添加数据，load文件		
	if(!optimizer.load("../data/added.g2o")){
		cout<<"Error loading graph"<<endl;
		return -1;
	}
	else{
		cout<<"Loaded "<<optimizer.vertices().size()<<" vertices"<<endl;
		cout<<"Loaded "<<optimizer.edges().size()<<" edges"<<endl;
	}
	
	//optimizer.edge->setRobustKernel( new RobustKernelHuber() );

	//优化过程中，第一个点固定，不做优化; 也可以不固定。
	VertexSE2* firstRobotPose = dynamic_cast<VertexSE2*>(optimizer.vertex(0));
	firstRobotPose->setFixed(true);	//这里固定了第一个点，如果固定最后一个点,或者不固定，结果会不一样
	
	optimizer.setAlgorithm(Algorithm);
	optimizer.setVerbose(true);
	optimizer.initializeOptimization();
	cout<<"Optimizing ..."<<endl;
	optimizer.optimize(MAXITERATION);

	cout<<"done."<<endl;

	optimizer.save("../data/optimized.g2o");
	optimizer.clear();

	return 0;  
}




void MyWrite(int num1, int num2, FILE *fpp){
	char str1[20], str3[20];
	sprintf(str1, "%d", num1);
	sprintf(str3, "%d", num2);
	char str0[200] = "EDGE_SE2 ";
	char str2[20] = " ";
	char str4[50] = " 0 0 0 900 0 0 1000 0 900\n";
	int cnt;
	
	strcat(str0, str1);
	strcat(str0, str2);
	strcat(str0, str3);
	strcat(str0, str4);
	cnt = fputs(str0, fpp);
	//printf("%s%d\n", str0, cnt);
}

float min3(float num1, float num2, float num3){
	if(num1 < num2){
		if(num1 < num3)
			return num1;
		else{
			if(num3 < num2)
				return num3;
		}
	}		
	else{
		if(num3<num2)
			return num3;
		else
			return num2;
	}
}

int min2(int num1, int num2){
	if(num1 < num2)
		return num1;
	else
		return num2;
}

void FileAnalysis(void){
	FILE *fp;
	char OneRaw[200];
	char temp;
	int read;
	char file[] = "../data/input.g2o";
	if((fp=fopen(file,"r+")) == NULL){
		perror("fail to open input.g2o!");
		exit(1);
	}
	else printf("open succeed input.g2o!\n");
	int point;
	fseek(fp, 0, SEEK_END);	//读写位置在文件末尾
	point = ftell(fp);
	printf("文件总字节数：%d\n", point);
	
	int last_cur=0;	//上一次的光标位置
	int iCnt;
	int row=0;	//记录总的行数
	char *delim = " ";
	char *p;
	for(iCnt=0;iCnt<point;iCnt++){	
		fseek(fp, iCnt, SEEK_SET);	//设定读取位置为i处
		temp = (char)fgetc(fp);	//读取一个字符		
		if(temp == '\n'){
			row++;
			fseek(fp, last_cur, SEEK_SET);	//设定读取位置为last_cur处	
			fgets(OneRaw, iCnt-last_cur, fp);	//
			last_cur = iCnt+1;
			//现在，OneRaw中存储的是一行的数据，接下来就可以对数据进行处理了
			if(OneRaw[0] == 'V'){	//找到了vertex的行	
				NumOfVertex++;
				printf("%s\n",OneRaw);
				//printf("%d\n",(int)sizeof(OneRaw));		
				p = strtok(OneRaw, delim);
				float num[4];
				int jCnt=0;
				while(p=strtok(NULL,delim)){
					num[jCnt] = atof(p);
					jCnt++;
					//printf("%s\n", p);
				}
				vertex_all[(unsigned int)num[0]].id = (unsigned int)num[0];
				vertex_all[(unsigned int)num[0]].x = num[1];
				vertex_all[(unsigned int)num[0]].y = num[2];
				vertex_all[(unsigned int)num[0]].theta = num[3];
				printf("%d\n",vertex_all[(unsigned int)num[0]].id);
				printf("%f\n",vertex_all[(unsigned int)num[0]].x);
				printf("%f\n",vertex_all[(unsigned int)num[0]].y);
				printf("%f\n",vertex_all[(unsigned int)num[0]].theta);
				//break;
			}
			else if(OneRaw[0] == 'E'){
				printf("%s\n",OneRaw);
				//printf("%d\n",(int)sizeof(OneRaw));		
				p = strtok(OneRaw, delim);
				float num[11];
				int jCnt=0;
				while(p=strtok(NULL,delim)){	//字符串分解
					num[jCnt] = atof(p);
					jCnt++;
					//printf("%s\n", p);
				}
				if((unsigned int)num[1]+1 == (unsigned int)num[0]){	//相邻边
					NumOfedge_Con++;
					edge_con[(unsigned int)num[1]].fromid = (unsigned int)num[0];
					edge_con[(unsigned int)num[1]].toid = (unsigned int)num[1];
					edge_con[(unsigned int)num[1]].delta_x = num[2];
					edge_con[(unsigned int)num[1]].delta_y = num[3];
					edge_con[(unsigned int)num[1]].delta_t = num[4];
					edge_con[(unsigned int)num[1]].info_xx = num[5];
					edge_con[(unsigned int)num[1]].info_xy = num[6];
					edge_con[(unsigned int)num[1]].info_xt = num[7];
					edge_con[(unsigned int)num[1]].info_yy = num[8];
					edge_con[(unsigned int)num[1]].info_yt = num[9];
					edge_con[(unsigned int)num[1]].info_tt = num[10];
					printf("%d\n",edge_con[(unsigned int)num[1]].fromid);
					printf("%d\n",edge_con[(unsigned int)num[1]].toid);
					printf("%f\n",edge_con[(unsigned int)num[1]].delta_x);
					printf("%f\n",edge_con[(unsigned int)num[1]].delta_y);
					printf("%f\n",edge_con[(unsigned int)num[1]].delta_t);
					printf("%f\n",edge_con[(unsigned int)num[1]].info_xx);
					printf("%f\n",edge_con[(unsigned int)num[1]].info_xy);
					printf("%f\n",edge_con[(unsigned int)num[1]].info_xt);
					printf("%f\n",edge_con[(unsigned int)num[1]].info_yy);
					printf("%f\n",edge_con[(unsigned int)num[1]].info_yt);
					printf("%f\n",edge_con[(unsigned int)num[1]].info_tt);
				}
				else{	//非相邻边
					NumOfedge_Ncon++;
					static int kCnt = 0;
					edge_Ncon[kCnt].fromid = (unsigned int)num[0];
					edge_Ncon[kCnt].toid = (unsigned int)num[1];
					edge_Ncon[kCnt].delta_x = num[2];
					edge_Ncon[kCnt].delta_y = num[3];
					edge_Ncon[kCnt].delta_t = num[4];
					edge_Ncon[kCnt].info_xx = num[5];
					edge_Ncon[kCnt].info_xy = num[6];
					edge_Ncon[kCnt].info_xt = num[7];
					edge_Ncon[kCnt].info_yy = num[8];
					edge_Ncon[kCnt].info_yt = num[9];
					edge_Ncon[kCnt].info_tt = num[10];
					printf("%d\n",kCnt);
					printf("%d\n",edge_Ncon[kCnt].fromid);
					printf("%d\n",edge_Ncon[kCnt].toid);
					printf("%f\n",edge_Ncon[kCnt].delta_x);
					printf("%f\n",edge_Ncon[kCnt].delta_y);
					printf("%f\n",edge_Ncon[kCnt].delta_t);
					printf("%f\n",edge_Ncon[kCnt].info_xx);
					printf("%f\n",edge_Ncon[kCnt].info_xy);
					printf("%f\n",edge_Ncon[kCnt].info_xt);
					printf("%f\n",edge_Ncon[kCnt].info_yy);
					printf("%f\n",edge_Ncon[kCnt].info_yt);
					printf("%f\n",edge_Ncon[kCnt].info_tt);
					kCnt++;
				}
				//break;
			}
		}
	}
	printf("节点数量：%d; 相邻边数量：%d; 非相邻边数量：%d;\n", NumOfVertex, NumOfedge_Con, NumOfedge_Ncon);
}

void MissCheck(void){	
	int qCnt;	
	//下面这段代码为数据分析，分析出了磁钉漏检的情况，可以作为程序设计时候的参考，结果在“磁钉漏检.txt中”	
	//分析方法：两圈测量的值差距在0.6m以上，则认为这个点测量错误。
	for(qCnt=0; qCnt<=min2(Cir1_Point2-Cir1_Point1, Cir2_Point2-Cir2_Point1); qCnt++){
		if(fabs(edge_con[Cir1_Point1+qCnt].delta_x - edge_con[Cir2_Point1+qCnt].delta_x)>0.6){
			printf("%d,  %d\n", Cir1_Point1+qCnt, Cir2_Point1+qCnt);
		}	
	}
	for(qCnt=0; qCnt<=min2(Cir1_Point3-Cir1_Point2, Cir2_Point3-Cir2_Point2); qCnt++){
		if(fabs(edge_con[Cir1_Point2+qCnt].delta_x - edge_con[Cir2_Point2+qCnt].delta_x)>0.6){
			printf("%d,  %d\n", Cir1_Point2+qCnt, Cir2_Point2+qCnt);
		}	
	}
	for(qCnt=0; qCnt<=min2(Cir1_Point4-Cir1_Point3, Cir2_Point4-Cir2_Point3); qCnt++){
		if(fabs(edge_con[Cir1_Point3+qCnt].delta_x - edge_con[Cir2_Point3+qCnt].delta_x)>0.6){
			printf("%d,  %d\n", Cir1_Point3+qCnt, Cir2_Point3+qCnt);
		}	
	}
	for(qCnt=0; qCnt<=min2(Cir1_Point5-Cir1_Point4, Cir2_Point5-Cir2_Point4); qCnt++){
		if(fabs(edge_con[Cir1_Point4+qCnt].delta_x - edge_con[Cir2_Point4+qCnt].delta_x)>0.6){
			printf("%d,  %d\n", Cir1_Point4+qCnt, Cir2_Point4+qCnt);
		}	
	}
//-----------------------------------------------------------------------------------------------------------	
	for(qCnt=0; qCnt<=min2(Cir2_Point2-Cir2_Point1, Cir3_Point2-Cir3_Point1); qCnt++){
		if(fabs(edge_con[Cir2_Point1+qCnt].delta_x - edge_con[Cir3_Point1+qCnt].delta_x)>0.6){
			printf("%d,  %d\n", Cir2_Point1+qCnt, Cir3_Point1+qCnt);
		}	
	}
	for(qCnt=0; qCnt<=min2(Cir2_Point3-Cir2_Point2, Cir3_Point3-Cir3_Point2); qCnt++){
		if(fabs(edge_con[Cir2_Point2+qCnt].delta_x - edge_con[Cir3_Point2+qCnt].delta_x)>0.6){
			printf("%d,  %d\n", Cir2_Point2+qCnt, Cir3_Point2+qCnt);
		}	
	}
	for(qCnt=0; qCnt<=min2(Cir2_Point4-Cir2_Point3, Cir3_Point4-Cir3_Point3); qCnt++){
		if(fabs(edge_con[Cir2_Point3+qCnt].delta_x - edge_con[Cir3_Point3+qCnt].delta_x)>0.6){
			printf("%d,  %d\n", Cir2_Point3+qCnt, Cir3_Point3+qCnt);
		}	
	}
	for(qCnt=0; qCnt<=min2(Cir2_Point5-Cir2_Point4, Cir3_Point5-Cir3_Point4); qCnt++){
		if(fabs(edge_con[Cir2_Point4+qCnt].delta_x - edge_con[Cir3_Point4+qCnt].delta_x)>0.6){
			printf("%d,  %d\n", Cir2_Point4+qCnt, Cir3_Point4+qCnt);
		}	
	}
	//上面这段代码为数据分析，分析出了磁钉漏检的情况，可以作为程序设计时候的参考，结果在“磁钉漏检.txt中”
}





