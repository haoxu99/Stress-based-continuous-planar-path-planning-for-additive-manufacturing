#include<iostream>
#include<OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>  
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <OpenMesh\Core\Utils\Property.hh> 
#include <string>
#include <vector>
#include <algorithm>
#include<stdio.h>
#include "clipper.hpp"  
#include"clipper.cpp"
#include<GL/glut.h>
//#include"GA.h"
#define use_int32
#define pi 3.1415926
using namespace std;
using namespace ClipperLib;

#define M 700
#define Die 5000
int  die = Die;

using namespace std;
struct group {
	vector<int>id;
	double p;
	double fit;
	double sum_p;
};
vector<vector<double>>dis;//距离矩阵
vector<int> groupbest;//每层的最优染色体
vector<group>groups;//染色体数组
double d1=0.56,d2=0.4;//路径宽度

int generation;
double groupbestp;//最优解的p
double groupbestfit;//最优解的fit
int changebest;//是否用最优解代替新种群
double distancenum;



struct element {//单元结构体，里面存储了每个单元的全部信息
	int id;
	double x, y, z, s11, s22, s12, s_max, direction;
	int b;
};
struct points {//存储bool为1 的单元信息
	int id;
	double x, y, z, s_max, direction;
	int b = 0;//该单元的度数
};
struct path {//路径结构体，每条路径都是一个双向
	int id;
	double x, y, z, direction;
	int b = 0;//判断该路径点是否为出入点，其中出点为-1，一般的点为0,入点为1
};



vector<element>point_buffer1;//第一层缓存
vector<vector<element>>point_buffer2;//第二层缓存
vector<vector<vector<element>>>point;//每个单元的信息
vector<vector<points>>model;//按应力大小，对有效单元进行排序，大->小,一维代表层数
vector<points>model_buffer;
vector<vector<vector<path>>>lines;
vector<vector<vector<path>>>paths;//存储路径，一维代表层数
vector<path>paths_buffer1;
vector<vector<path>>paths_buffer2;
vector<vector<vector<path>>>paths_buffer3;
vector<vector<vector<path>>>path_w;//外轮廓
queue<points> que;

vector<vector<vector<int>>>G;//邻接矩阵，一维为层
vector<int>G_buffer1;
vector<vector<int>>G_buffer2;
vector<vector<int>>visited;//访问数组
vector<int>arr_buffer1;
vector<vector<int>>arr;//存储了每层路径的遍历循序
//int num=0;//工具
int x_count, y_count, z_count;//分别表示x，y,z的单元个数
double Bmin_x, Bmin_y, Bmax_x, Bmax_y;

void ReadFile() {
	element e;
	ifstream file1, file2,file3;
	file1.open("VoxelizationData.txt");//读取体素文件
	file2.open("abaqus.txt");//读取应力文件
	file3.open("path_w.txt");
	double a1, a2;
	int q1 = 25095;//行数
	path s;
	/*for (int i = 0; i < q1; i++) {
		file3 >> a1;
		file3 >> a2;
		if (a1 != 9999 && a2 != 9999) {
			s.x = a1;
			s.y = a2;
			paths_buffer1.push_back(s);
		}
		else {
			path_w.push_back(paths_buffer1);
			paths_buffer1.clear();
		}
	}
	path_w.push_back(paths_buffer1);
	paths_buffer1.clear();*/
	for (int i = 0; i < q1; i++) {
		file3 >> a1;
		file3 >> a2;
		if ((a1 != 8888 && a2 != 8888) && (a1 != 9999 && a2 != 9999)) {
			s.x = a1;
			s.y = a2;
			paths_buffer1.push_back(s);
		}
		else if (a1 == 9999 && a2 == 9999) {
			paths_buffer2.push_back(paths_buffer1);
			paths_buffer1.clear();
		}
		else {
			paths_buffer2.push_back(paths_buffer1);
			paths_buffer1.clear();
			path_w.push_back(paths_buffer2);
			paths_buffer2.clear();
		}
	}
	paths_buffer1.clear();
	paths_buffer2.clear();

	file1 >> x_count;
	file1 >> y_count;
	file1 >> z_count;
	printf("%d %d %d \n", x_count, y_count, z_count);
	for (int i = 0; i < z_count; i++) {//将每个单元的坐标信息存入数组
		for (int j = 0; j < y_count; j++) {
			for (int m = 0; m < x_count; m++) {
				for (int s = 0; s < 5; s++) {
					switch (s) {
					case 0:file1 >> e.id; break;
					case 1:file1 >> e.x; break;
					case 2:file1 >> e.y; break;
					case 3:file1 >> e.z; break;
					case 4:file1 >> e.b; break;
					}
				}
				point_buffer1.push_back(e);
			}
			point_buffer2.push_back(point_buffer1);
			point_buffer1.clear();
		}
		point.push_back(point_buffer2);
		point_buffer2.clear();
	}
	for (int i = 0; i < point.size(); i++) {//读取单元的x,y应力数据
		for (int j = 0; j < point[0].size(); j++) {
			for (int m = 0; m < point[0][0].size(); m++) {
				if (point[i][j][m].b == 1) {
					for (int s = 0; s < 4; s++) {
						switch (s) {
						case 0: {file2 >> e.id; if (e.id != point[i][j][m].id) printf("文件读取错误！\n"); break; }
						case 1:file2 >> point[i][j][m].s11; break;
						case 2:file2 >> point[i][j][m].s22; break;
						case 3:file2 >> point[i][j][m].s12; break;
						}
					}
				}
			}
		}
	}
}

double distance(path a, path b)	//两点间距离
{
	return sqrt(pow((a.x - b.x), 2) + pow((a.y - b.y), 2));
}

void FE_Analysis() {//计算最大主应力以及最大主应力方向
	double s1, s2, d1, d2;
	for (int i = 0; i < point.size(); i++) {
		for (int j = 0; j < point[i].size(); j++) {
			for (int m = 0; m < point[i][j].size(); m++) {
				if (point[i][j][m].b == 1) {
					//最大主应力数值（abs）
					s1 = ((point[i][j][m].s11 + point[i][j][m].s22) / (double)2) - sqrt(pow(((point[i][j][m].s11 - point[i][j][m].s22) / (double)2), 2) + pow(point[i][j][m].s12, 2));
					s2 = ((point[i][j][m].s11 + point[i][j][m].s22) / (double)2) + sqrt(pow(((point[i][j][m].s11 - point[i][j][m].s22) / (double)2), 2) + pow(point[i][j][m].s12, 2));
					if (fabs(s1) >= fabs(s2)) {
						point[i][j][m].s_max = s1;
					}
					else {
						point[i][j][m].s_max = s2;
					}
					//最大主应力方向
					d1 = atan(-(2 * point[i][j][m].s12) / (point[i][j][m].s11 - point[i][j][m].s22));
					if (point[i][j][m].s11 >= point[i][j][m].s22) {
						if (d1 <= 0) d1 = d1 + 2 * pi;
						
					}
					else if (point[i][j][m].s11 <= point[i][j][m].s22) {
						d1 = d1 + pi;
					}
					d1 = d1 / (double)2;
					//最小主应力方向
					d2= atan(-(2 * point[i][j][m].s12) / (point[i][j][m].s11 - point[i][j][m].s22));
					if (point[i][j][m].s11 >= point[i][j][m].s22) {
						d2 = d2 + pi;
					}
					else if (point[i][j][m].s11 <= point[i][j][m].s11) {
						if (d2 <= 0) d2 = d2 + 2 * pi;
					}
					d2 = d2 / (double)2;
					//最大主应力（abs）方向
					if ((s1 > s2 && fabs(s1) > fabs(s2)) || (s2 > s1 && fabs(s2) > fabs(s1))) point[i][j][m].direction = pi-d1;
					else if ((s1 > s2 && fabs(s1) < fabs(s2)) || (s2 > s1 && fabs(s2) < fabs(s1))) point[i][j][m].direction =pi- d2;
				}
			}
		}
	}
}

void VectorInsert(int i,int m) {//生成单路径函数
	bool flag = true;
	path s;
	int Right_up, Left_up, Right_down, Left_down;
	int q=0;
	while (flag) {
		q++;
		flag = false;
		Right_up = 0;
		Left_up = 0;
		Right_down = 0;
		Left_down = 0;
			for (int c = 0; c < paths_buffer2.size(); c++) {
					for (int n = 1; n < paths_buffer2[c].size(); n++) {
						if ((model[i][m].id + 1 == paths_buffer2[c][n - 1].id && model[i][m].id + x_count == paths_buffer2[c][n].id && model[i][m].id % x_count != 0) ||
							(model[i][m].id + 1 == paths_buffer2[c][n].id && model[i][m].id + x_count == paths_buffer2[c][n - 1].id && model[i][m].id % x_count != 0)) {
							Right_up = 1;//如果与右上的端点连接，则会与已有的路径交叉
						}
						if ((model[i][m].id - 1 == paths_buffer2[c][n - 1].id && model[i][m].id + x_count == paths_buffer2[c][n].id && model[i][m].id % x_count != 1) ||
							(model[i][m].id - 1 == paths_buffer2[c][n].id && model[i][m].id + x_count == paths_buffer2[c][n - 1].id && model[i][m].id % x_count != 1)) {
							Left_up = 1;//如果与左上的端点连接，则会与已有的路径交叉
						}
						if ((model[i][m].id + 1 == paths_buffer2[c][n - 1].id && model[i][m].id - x_count == paths_buffer2[c][n].id && model[i][m].id % x_count != 0) ||
							(model[i][m].id + 1 == paths_buffer2[c][n].id && model[i][m].id - x_count == paths_buffer2[c][n - 1].id && model[i][m].id % x_count != 0)) {
							Right_down = 1;//如果与右下的端点连接，则会与已有的路径交叉
						}
						if ((model[i][m].id - 1 == paths_buffer2[c][n - 1].id && model[i][m].id - x_count == paths_buffer2[c][n].id && model[i][m].id % x_count != 1) ||
							(model[i][m].id - 1 == paths_buffer2[c][n].id && model[i][m].id - x_count == paths_buffer2[c][n - 1].id && model[i][m].id % x_count != 1)) {
							Left_down = 1;//如果与左下的端点连接，则会与已有的路径交叉
						}
					}
			}
		if (model[i][m].direction >= 0 && model[i][m].direction < pi / (double)8) {
			
			for (int j = 0; j < model[i].size(); j++) {
				if (((model[i][j].id == model[i][m].id + 1 && model[i][m].id % x_count != 0) || (model[i][j].id == model[i][m].id - 1 && model[i][m].id % x_count != 1)) && model[i][j].b == 0) {//右边的单元
					model[i][m].b++;
					model[i][j].b++;
					m = j;
					s.id = model[i][j].id;
					s.x = model[i][j].x;
					s.y = model[i][j].y;
					s.z = model[i][j].z;
					s.direction = model[i][j].direction;
					paths_buffer1.push_back(s);
					//flag = true;
					if (model[i][m].direction >= 0 && model[i][m].direction < pi / (double)8)
					{
						flag = true;
					}
					break;
				}
			}
		}
		else if (model[i][m].direction >= pi / (double)8 && model[i][m].direction < pi * ((double)3 / (double)8)) {
			for (int j = 0; j < model[i].size(); j++) {
				if (((model[i][j].id == model[i][m].id + x_count + 1 && model[i][m].id % x_count != 0 && Right_up == 0) || (model[i][j].id == model[i][m].id - x_count - 1 && model[i][m].id % x_count != 1 && Left_down == 0))  && model[i][j].b == 0) {//右上边的单元
					model[i][m].b++;
					model[i][j].b++;
					m = j;
					s.id = model[i][j].id;
					s.x = model[i][j].x;
					s.y = model[i][j].y;
					s.z = model[i][j].z;
					s.direction = model[i][j].direction;
					paths_buffer1.push_back(s);
					//flag = true;
					if (model[i][m].direction >= pi / (double)8 && model[i][m].direction < pi * ((double)3 / (double)8))
					{
						flag = true;
					}
					
					break;
				}
			}
		}
		else if (model[i][m].direction >= pi * ((double)3 / (double)8) && model[i][m].direction < pi * ((double)5 / (double)8)) {
			for (int j = 0; j < model[i].size(); j++) {
				if ((model[i][j].id == model[i][m].id + x_count || model[i][j].id == model[i][m].id - x_count) && model[i][j].b == 0) {//上边的单元
					model[i][m].b++;
					model[i][j].b++;
					m = j;
					s.id = model[i][j].id;
					s.x = model[i][j].x;
					s.y = model[i][j].y;
					s.z = model[i][j].z;
					s.direction = model[i][j].direction;
					paths_buffer1.push_back(s);
					//flag = true;
					if (model[i][m].direction >= pi * ((double)3 / (double)8) && model[i][m].direction < pi * ((double)5 / (double)8))
					{
						flag = true;
					}
					
					break;
				}
			}
		}
		else if (model[i][m].direction >= pi * ((double)5 / (double)8) && model[i][m].direction < pi * ((double)7 / (double)8)) {
			for (int j = 0; j < model[i].size(); j++) {
				if (((model[i][j].id == model[i][m].id + x_count - 1 && model[i][m].id % x_count != 1 && Left_up == 0) || (model[i][j].id == model[i][m].id - x_count + 1 && model[i][m].id % x_count != 0 && Right_down == 0)) && model[i][j].b == 0) {//左上边的单元
					model[i][m].b++;
					model[i][j].b++;
					m = j;
					s.id = model[i][j].id;
					s.x = model[i][j].x;
					s.y = model[i][j].y;
					s.z = model[i][j].z;
					s.direction = model[i][j].direction;
					paths_buffer1.push_back(s);
					//flag = true;
					if (model[i][m].direction >= pi * ((double)5 / (double)8) && model[i][m].direction < pi * ((double)7 / (double)8))
					{
						flag = true;
					}
					
					break;
				}
			}
		}
		else if (model[i][m].direction >= pi * ((double)7 / (double)8) && model[i][m].direction < pi) {
			for (int j = 0; j < model[i].size(); j++) {
				if (((model[i][j].id == model[i][m].id - 1 && model[i][m].id % x_count != 1) || (model[i][j].id == model[i][m].id + 1 && model[i][m].id % x_count != 0)) && model[i][j].b == 0) {//左边的单元
					model[i][m].b++;
					model[i][j].b++;
					m = j;
					s.id = model[i][j].id;
					s.x = model[i][j].x;
					s.y = model[i][j].y;
					s.z = model[i][j].z;
					s.direction = model[i][j].direction;
					paths_buffer1.push_back(s);
					//flag = true;
					if (model[i][m].direction >= pi * ((double)7 / (double)8) && model[i][m].direction < pi)
					{
						flag = true;
					}
					
					break;
				}
			}
		}
	}
}

void OptimizePoints() {//优化空点函数
	path s;
	int Right_up = 0;//判断该空点的斜连接是否可连接，可以连接为0，不可连接为1
	int Left_up = 0;
	int Right_down = 0;
	int Left_down = 0;
	for (int i = 0; i < model.size(); i++) {
		for (int j = 0; j < model[i].size(); j++) {
			if (model[i][j].b == 0) {
				Right_up = 0;//判断该空点的斜连接是否可连接，可以连接为0，不可连接为1
				Left_up = 0;
				Right_down = 0;
				Left_down = 0;
				for (int m = 0; m < paths[i].size(); m++) {
					for (int n = 1; n < paths[i][m].size(); n++) {
						if ((model[i][j].id + 1 == paths[i][m][n - 1].id && model[i][j].id + x_count == paths[i][m][n].id && model[i][j].id % x_count != 0) ||
							(model[i][j].id + 1 == paths[i][m][n].id && model[i][j].id + x_count == paths[i][m][n - 1].id && model[i][j].id % x_count != 0)) {
							Right_up = 1;//如果与右上的端点连接，则会与已有的路径交叉
						}
						if ((model[i][j].id - 1 == paths[i][m][n - 1].id && model[i][j].id + x_count == paths[i][m][n].id && model[i][j].id % x_count != 1) ||
							(model[i][j].id - 1 == paths[i][m][n].id && model[i][j].id + x_count == paths[i][m][n - 1].id && model[i][j].id % x_count != 1)) {
							Left_up = 1;//如果与左上的端点连接，则会与已有的路径交叉
						}
						if ((model[i][j].id + 1 == paths[i][m][n - 1].id && model[i][j].id - x_count == paths[i][m][n].id && model[i][j].id % x_count != 0) ||
							(model[i][j].id + 1 == paths[i][m][n].id && model[i][j].id - x_count == paths[i][m][n - 1].id && model[i][j].id % x_count != 0)) {
							Right_down = 1;//如果与右下的端点连接，则会与已有的路径交叉
						}
						if ((model[i][j].id - 1 == paths[i][m][n - 1].id && model[i][j].id - x_count == paths[i][m][n].id && model[i][j].id % x_count != 1) ||
							(model[i][j].id - 1 == paths[i][m][n].id && model[i][j].id - x_count == paths[i][m][n - 1].id && model[i][j].id % x_count != 1)) {
							Left_down = 1;//如果与左下的端点连接，则会与已有的路径交叉
						}
					}
				}

				for (int m = 0; m < paths[i].size(); m++) {
					//右-终点
					if (paths[i][m][paths[i][m].size() - 2].id == paths[i][m][paths[i][m].size() - 1].id + 1 && paths[i][m][paths[i][m].size() - 1].id == model[i][j].id + 1 && model[i][j].id % x_count != 0) {
						model[i][j].b++;
						s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
						paths[i][m].push_back(s);
						j = 0;
						break;//寻找到连接的路径，退出该空点的循环
					}//上-终点
					else if (paths[i][m][paths[i][m].size() - 2].id == paths[i][m][paths[i][m].size() - 1].id + x_count && paths[i][m][paths[i][m].size() - 1].id == model[i][j].id + x_count) {
						model[i][j].b++;
						s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
						paths[i][m].push_back(s);
						j = 0;
						break;
					}//左-终点
					else if (paths[i][m][paths[i][m].size() - 2].id == paths[i][m][paths[i][m].size() - 1].id - 1 && paths[i][m][paths[i][m].size() - 1].id == model[i][j].id - 1 && model[i][j].id % x_count != 1) {
						model[i][j].b++;
						s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
						paths[i][m].push_back(s);
						j = 0;
						break;
					}//下-终点
					else if (paths[i][m][paths[i][m].size() - 2].id == paths[i][m][paths[i][m].size() - 1].id - x_count && paths[i][m][paths[i][m].size() - 1].id == model[i][j].id - x_count) {
						model[i][j].b++;
						s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
						paths[i][m].push_back(s);
						j = 0;
						break;
					}//右上-终点
					else if (paths[i][m][paths[i][m].size() - 2].id == paths[i][m][paths[i][m].size() - 1].id + 1 + x_count && paths[i][m][paths[i][m].size() - 1].id == model[i][j].id + 1 + x_count && model[i][j].id % x_count != 0 && Right_up == 0) {
						model[i][j].b++;
						s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
						paths[i][m].push_back(s);
						j = 0;
						break;
					}//左上-终点
					else if (paths[i][m][paths[i][m].size() - 2].id == paths[i][m][paths[i][m].size() - 1].id - 1 + x_count && paths[i][m][paths[i][m].size() - 1].id == model[i][j].id - 1 + x_count && model[i][j].id % x_count != 1 && Left_up == 0) {
						model[i][j].b++;
						s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
						paths[i][m].push_back(s);
						j = 0;
						break;
					}//左下-终点
					else if (paths[i][m][paths[i][m].size() - 2].id == paths[i][m][paths[i][m].size() - 1].id - 1 - x_count && paths[i][m][paths[i][m].size() - 1].id == model[i][j].id - 1 - x_count && model[i][j].id % x_count != 1 && Left_down == 0) {
						model[i][j].b++;
						s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
						paths[i][m].push_back(s);
						j = 0;
						break;
					}//右下-终点
					else if (paths[i][m][paths[i][m].size() - 2].id == paths[i][m][paths[i][m].size() - 1].id + 1 - x_count && paths[i][m][paths[i][m].size() - 1].id == model[i][j].id + 1 - x_count && model[i][j].id % x_count != 0 && Right_down == 0) {
						model[i][j].b++;
						s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
						paths[i][m].push_back(s);
						j = 0;
						break;
					}//右-起点
					else if (paths[i][m][1].id == paths[i][m][0].id + 1 && paths[i][m][0].id == model[i][j].id + 1 && model[i][j].id % x_count != 0) {
						model[i][j].b++;
						s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
						paths[i][m].push_back(s);
						for (int f = paths[i][m].size() - 1; f > 0; f--) {
							s = paths[i][m][f];
							paths[i][m][f] = paths[i][m][f - 1];
							paths[i][m][f - 1] = s;
						}
						j = 0;
						break;
					}//上-起点
					else if (paths[i][m][1].id == paths[i][m][0].id + x_count && paths[i][m][0].id == model[i][j].id + x_count) {
						model[i][j].b++;
						s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
						paths[i][m].push_back(s);
						for (int f = paths[i][m].size() - 1; f > 0; f--) {
							s = paths[i][m][f];
							paths[i][m][f] = paths[i][m][f - 1];
							paths[i][m][f - 1] = s;
						}
						j = 0;
						break;
					}//左-起点
					else if (paths[i][m][1].id == paths[i][m][0].id - 1 && paths[i][m][0].id == model[i][j].id - 1 && model[i][j].id % x_count != 1) {
						model[i][j].b++;
						s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
						paths[i][m].push_back(s);
						for (int f = paths[i][m].size() - 1; f > 0; f--) {
							s = paths[i][m][f];
							paths[i][m][f] = paths[i][m][f - 1];
							paths[i][m][f - 1] = s;
						}
						j = 0;
						break;
					}//下-起点
					else if (paths[i][m][1].id == paths[i][m][0].id - x_count && paths[i][m][0].id == model[i][j].id - x_count) {
						model[i][j].b++;
						s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
						paths[i][m].push_back(s);
						for (int f = paths[i][m].size() - 1; f > 0; f--) {
							s = paths[i][m][f];
							paths[i][m][f] = paths[i][m][f - 1];
							paths[i][m][f - 1] = s;
						}
						j = 0;
						break;
					}//右上-起点
					else if (paths[i][m][1].id == paths[i][m][0].id + 1 + x_count && paths[i][m][0].id == model[i][j].id + 1 + x_count && model[i][j].id % x_count != 0 && Right_up == 0) {
						model[i][j].b++;
						s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
						paths[i][m].push_back(s);
						for (int f = paths[i][m].size() - 1; f > 0; f--) {
							s = paths[i][m][f];
							paths[i][m][f] = paths[i][m][f - 1];
							paths[i][m][f - 1] = s;
						}
						j = 0;
						break;
					}//左上-起点
					else if (paths[i][m][1].id == paths[i][m][0].id - 1 + x_count && paths[i][m][0].id == model[i][j].id - 1 + x_count && model[i][j].id % x_count != 1 && Left_up == 0) {
						model[i][j].b++;
						s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
						paths[i][m].push_back(s);
						for (int f = paths[i][m].size() - 1; f > 0; f--) {
							s = paths[i][m][f];
							paths[i][m][f] = paths[i][m][f - 1];
							paths[i][m][f - 1] = s;
						}
						j = 0;
						break;
					}//左下-起点
					else if (paths[i][m][1].id == paths[i][m][0].id - 1 - x_count && paths[i][m][0].id == model[i][j].id - 1 - x_count && model[i][j].id % x_count != 1 && Left_down == 0) {
						model[i][j].b++;
						s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
						paths[i][m].push_back(s);
						for (int f = paths[i][m].size() - 1; f > 0; f--) {
							s = paths[i][m][f];
							paths[i][m][f] = paths[i][m][f - 1];
							paths[i][m][f - 1] = s;
						}
						j = 0;
						break;
					}//右下-起点
					else if (paths[i][m][1].id == paths[i][m][0].id + 1 - x_count && paths[i][m][0].id == model[i][j].id + 1 - x_count && model[i][j].id % x_count != 0 && Right_down == 0) {
						model[i][j].b++;
						s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
						paths[i][m].push_back(s);
						for (int f = paths[i][m].size() - 1; f > 0; f--) {
							s = paths[i][m][f];
							paths[i][m][f] = paths[i][m][f - 1];
							paths[i][m][f - 1] = s;
						}
						j = 0;
						break;
					}


					////右-终点
					//if (paths[i][m][paths[i][m].size() - 1].id == model[i][j].id + 1 && model[i][j].id % x_count != 0) {
					//	model[i][j].b++;
					//	s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					//	paths[i][m].push_back(s);
					//	j = 0;
					//	break;//寻找到连接的路径，退出该空点的循环
					//}//上-终点
					//else if (paths[i][m][paths[i][m].size() - 1].id == model[i][j].id + x_count) {
					//	model[i][j].b++;
					//	s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					//	paths[i][m].push_back(s);
					//	j = 0;
					//	break;
					//}//左-终点
					//else if (paths[i][m][paths[i][m].size() - 1].id == model[i][j].id - 1 && model[i][j].id % x_count != 1) {
					//	model[i][j].b++;
					//	s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					//	paths[i][m].push_back(s);
					//	j = 0;
					//	break;
					//}//下-终点
					//else if (paths[i][m][paths[i][m].size() - 1].id == model[i][j].id - x_count) {
					//	model[i][j].b++;
					//	s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					//	paths[i][m].push_back(s);
					//	j = 0;
					//	break;
					//}//右上-终点
					//else if (paths[i][m][paths[i][m].size() - 1].id == model[i][j].id + 1 + x_count && model[i][j].id % x_count != 0 && Right_up == 0) {
					//	model[i][j].b++;
					//	s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					//	paths[i][m].push_back(s);
					//	j = 0;
					//	break;
					//}//左上-终点
					//else if (paths[i][m][paths[i][m].size() - 1].id == model[i][j].id - 1 + x_count && model[i][j].id % x_count != 1 && Left_up == 0) {
					//	model[i][j].b++;
					//	s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					//	paths[i][m].push_back(s);
					//	j = 0;
					//	break;
					//}//左下-终点
					//else if (paths[i][m][paths[i][m].size() - 1].id == model[i][j].id - 1 - x_count && model[i][j].id % x_count != 1 && Left_down == 0) {
					//	model[i][j].b++;
					//	s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					//	paths[i][m].push_back(s);
					//	j = 0;
					//	break;
					//}//右下-终点
					//else if (paths[i][m][paths[i][m].size() - 1].id == model[i][j].id + 1 - x_count && model[i][j].id % x_count != 0 && Right_down == 0) {
					//	model[i][j].b++;
					//	s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					//	paths[i][m].push_back(s);
					//	j = 0;
					//	break;
					//}//右-起点
					//else if (paths[i][m][0].id == model[i][j].id + 1 && model[i][j].id % x_count != 0) {
					//	model[i][j].b++;
					//	s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					//	paths[i][m].push_back(s);
					//	for (int f = paths[i][m].size() - 1; f > 0; f--) {
					//		s = paths[i][m][f];
					//		paths[i][m][f] = paths[i][m][f - 1];
					//		paths[i][m][f - 1] = s;
					//	}
					//	j = 0;
					//	break;
					//}//上-起点
					//else if (paths[i][m][0].id == model[i][j].id + x_count) {
					//	model[i][j].b++;
					//	s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					//	paths[i][m].push_back(s);
					//	for (int f = paths[i][m].size() - 1; f > 0; f--) {
					//		s = paths[i][m][f];
					//		paths[i][m][f] = paths[i][m][f - 1];
					//		paths[i][m][f - 1] = s;
					//	}
					//	j = 0;
					//	break;
					//}//左-起点
					//else if (paths[i][m][0].id == model[i][j].id - 1 && model[i][j].id % x_count != 1) {
					//	model[i][j].b++;
					//	s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					//	paths[i][m].push_back(s);
					//	for (int f = paths[i][m].size() - 1; f > 0; f--) {
					//		s = paths[i][m][f];
					//		paths[i][m][f] = paths[i][m][f - 1];
					//		paths[i][m][f - 1] = s;
					//	}
					//	j = 0;
					//	break;
					//}//下-起点
					//else if (paths[i][m][0].id == model[i][j].id - x_count) {
					//	model[i][j].b++;
					//	s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					//	paths[i][m].push_back(s);
					//	for (int f = paths[i][m].size() - 1; f > 0; f--) {
					//		s = paths[i][m][f];
					//		paths[i][m][f] = paths[i][m][f - 1];
					//		paths[i][m][f - 1] = s;
					//	}
					//	j = 0;
					//	break;
					//}//右上-起点
					//else if (paths[i][m][0].id == model[i][j].id + 1 + x_count && model[i][j].id % x_count != 0 && Right_up == 0) {
					//	model[i][j].b++;
					//	s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					//	paths[i][m].push_back(s);
					//	for (int f = paths[i][m].size() - 1; f > 0; f--) {
					//		s = paths[i][m][f];
					//		paths[i][m][f] = paths[i][m][f - 1];
					//		paths[i][m][f - 1] = s;
					//	}
					//	j = 0;
					//	break;
					//}//左上-起点
					//else if (paths[i][m][0].id == model[i][j].id - 1 + x_count && model[i][j].id % x_count != 1 && Left_up == 0) {
					//	model[i][j].b++;
					//	s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					//	paths[i][m].push_back(s);
					//	for (int f = paths[i][m].size() - 1; f > 0; f--) {
					//		s = paths[i][m][f];
					//		paths[i][m][f] = paths[i][m][f - 1];
					//		paths[i][m][f - 1] = s;
					//	}
					//	j = 0;
					//	break;
					//}//左下-起点
					//else if (paths[i][m][0].id == model[i][j].id - 1 - x_count && model[i][j].id % x_count != 1 && Left_down == 0) {
					//	model[i][j].b++;
					//	s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					//	paths[i][m].push_back(s);
					//	for (int f = paths[i][m].size() - 1; f > 0; f--) {
					//		s = paths[i][m][f];
					//		paths[i][m][f] = paths[i][m][f - 1];
					//		paths[i][m][f - 1] = s;
					//	}
					//	j = 0;
					//	break;
					//}//右下-起点
					//else if (paths[i][m][0].id == model[i][j].id + 1 - x_count && model[i][j].id % x_count != 0 && Right_down == 0) {
					//	model[i][j].b++;
					//	s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					//	paths[i][m].push_back(s);
					//	for (int f = paths[i][m].size() - 1; f > 0; f--) {
					//		s = paths[i][m][f];
					//		paths[i][m][f] = paths[i][m][f - 1];
					//		paths[i][m][f - 1] = s;
					//	}
					//	j = 0;
					//	break;
					//}


				}
			}
		}
	}
}

void PathConnect() {
	paths_buffer1.clear();
	path s;
	int Right_up, Left_up, Right_down, Left_down;
	///////////////////////根据方向先进行合并（终点-起点）
	for (int i = 0; i < paths.size(); i++) {

	
	///////////同斜率的合并
	

		for (int j = 0; j < paths[i].size(); j++) {
			Right_up = 0;
			Left_up = 0;
			Right_down = 0;
			Left_down = 0;
			for (int c = 0; c < paths[i].size(); c++) {
				for (int n = 1; n < paths[i][c].size(); n++) {
					if ((paths[i][j][0].id + 1 == paths[i][c][n - 1].id && paths[i][j][0].id + x_count == paths[i][c][n].id && paths[i][j][0].id % x_count != 0) ||
						(paths[i][j][0].id + 1 == paths[i][c][n].id && paths[i][j][0].id + x_count == paths[i][c][n - 1].id && paths[i][j][0].id % x_count != 0)) {
						Right_up = 1;//如果与右上的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][0].id - 1 == paths[i][c][n - 1].id && paths[i][j][0].id + x_count == paths[i][c][n].id && paths[i][j][0].id % x_count != 1) ||
						(paths[i][j][0].id - 1 == paths[i][c][n].id && paths[i][j][0].id + x_count == paths[i][c][n - 1].id && paths[i][j][0].id % x_count != 1)) {
						Left_up = 1;//如果与左上的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][0].id + 1 == paths[i][c][n - 1].id && paths[i][j][0].id - x_count == paths[i][c][n].id && paths[i][j][0].id % x_count != 0) ||
						(paths[i][j][0].id + 1 == paths[i][c][n].id && paths[i][j][0].id - x_count == paths[i][c][n - 1].id && paths[i][j][0].id % x_count != 0)) {
						Right_down = 1;//如果与右下的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][0].id - 1 == paths[i][c][n - 1].id && paths[i][j][0].id - x_count == paths[i][c][n].id && paths[i][j][0].id % x_count != 1) ||
						(paths[i][j][0].id - 1 == paths[i][c][n].id && paths[i][j][0].id - x_count == paths[i][c][n - 1].id && paths[i][j][0].id % x_count != 1)) {
						Left_down = 1;//如果与左下的端点连接，则会与已有的路径交叉
					}
				}
			}
			for (int m = 0; m < paths[i].size(); m++) {

				if (m == j) m++;
				if (m == paths[i].size()) break;
				//起点-起点
				if (paths[i][m][1].id  == 2*paths[i][m][0].id - paths[i][j][0].id && paths[i][m][0].id ==  2*paths[i][j][0].id - paths[i][j][1].id) {
					if ((paths[i][m][0].id == paths[i][j][0].id + 1 + x_count && Right_up == 1) || (paths[i][m][0].id == paths[i][j][0].id - 1 + x_count && Left_up == 1) || (paths[i][m][0].id == paths[i][j][0].id - 1 - x_count && Left_down == 1) || (paths[i][m][0].id == paths[i][j][0].id + 1 - x_count && Right_down == 1))
						break;
					if ((paths[i][j][0].id % x_count == 0 && (paths[i][m][0].id == paths[i][j][0].id + 1 || paths[i][m][0].id == paths[i][j][0].id + 1 + x_count || paths[i][m][0].id == paths[i][j][0].id + 1 - x_count)) || (paths[i][j][0].id % x_count == 1 && (paths[i][m][0].id == paths[i][j][0].id - 1 || paths[i][m][0].id == paths[i][j][0].id - 1 + x_count || paths[i][m][0].id == paths[i][j][0].id - 1 - x_count)))
						break;
					reverse(paths[i][j].begin(), paths[i][j].end());//翻转
					paths[i][j].insert(paths[i][j].end(), paths[i][m].begin(), paths[i][m].end());
					paths[i].erase(paths[i].begin() + m);
					m = 0; j = 0;
					break;
				}
				//起点-终点
				if (paths[i][m][paths[i][m].size() - 2].id  == 2*paths[i][m][paths[i][m].size() - 1].id - paths[i][j][0].id && paths[i][m][paths[i][m].size() - 1].id == 2*paths[i][j][0].id - paths[i][j][1].id ) {
					if ((paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id + 1 + x_count && Right_up == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id - 1 + x_count && Left_up == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id - 1 - x_count && Left_down == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id + 1 - x_count && Right_down == 1))
						break;
					if ((paths[i][j][0].id % x_count == 0 && (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id + 1 || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id + 1 + x_count || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id + 1 - x_count)) || (paths[i][j][0].id % x_count == 1 && (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id - 1 || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id - 1 + x_count || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id - 1 - x_count)))
						break;
					paths[i][m].insert(paths[i][m].end(), paths[i][j].begin(), paths[i][j].end());
					paths[i].erase(paths[i].begin() + j);
					m = 0; j = 0;
					break;
				}
			}
		}

		

		for (int j = 0; j < paths[i].size(); j++) {
			Right_up = 0;//判断该空点的斜连接是否可连接
			//可以连接为0，不可连接为1
			Left_up = 0;
			Right_down = 0;
			Left_down = 0;
			for (int c = 0; c < paths[i].size(); c++) {
				for (int n = 1; n < paths[i][c].size(); n++) {
					if ((paths[i][j][paths[i][j].size() - 1].id + 1 == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id + x_count == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 0) ||
						(paths[i][j][paths[i][j].size() - 1].id + 1 == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id + x_count == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 0)) {
						Right_up = 1;//如果与右上的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][paths[i][j].size() - 1].id - 1 == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id + x_count == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 1) ||
						(paths[i][j][paths[i][j].size() - 1].id - 1 == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id + x_count == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 1)) {
						Left_up = 1;//如果与左上的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][paths[i][j].size() - 1].id + 1 == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id - x_count == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 0) ||
						(paths[i][j][paths[i][j].size() - 1].id + 1 == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id - x_count == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 0)) {
						Right_down = 1;//如果与右下的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][paths[i][j].size() - 1].id - 1 == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id - x_count == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 1) ||
						(paths[i][j][paths[i][j].size() - 1].id - 1 == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id - x_count == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 1)) {
						Left_down = 1;//如果与左下的端点连接，则会与已有的路径交叉
					}
				}
			}
			for (int m = 0; m < paths[i].size(); m++) {



				if (m == j) m++;
				if (m == paths[i].size()) break;
				////终点-终点
				if (paths[i][m][paths[i][m].size() - 2].id==2*paths[i][m][paths[i][m].size() - 1].id - paths[i][j][paths[i][j].size() - 1].id && paths[i][m][paths[i][m].size() - 1].id== 2*paths[i][j][paths[i][j].size() - 1].id - paths[i][j][paths[i][j].size() - 2].id ) {
					if ((paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id + 1 + x_count && Right_up == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id - 1 + x_count && Left_up == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id - 1 - x_count && Left_down == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id + 1 - x_count && Right_down == 1))
						break;
					if ((paths[i][j][paths[i][j].size() - 1].id % x_count == 0 && (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id + 1 || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id + 1 + x_count || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id + 1 - x_count)) || (paths[i][j][paths[i][j].size() - 1].id % x_count == 1 && (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id - 1 || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id - 1 + x_count || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id - 1 - x_count)))
						break;
					reverse(paths[i][m].begin(), paths[i][m].end());//翻转
					paths[i][j].insert(paths[i][j].end(), paths[i][m].begin(), paths[i][m].end());
					paths[i].erase(paths[i].begin() + m);
					j = 0; m = 0;
					break;
				}
				////终点-起点
				if (paths[i][m][1].id ==2* paths[i][m][0].id - paths[i][j][paths[i][j].size() - 1].id && paths[i][m][0].id== 2*paths[i][j][paths[i][j].size() - 1].id - paths[i][j][paths[i][j].size() - 2].id ) {
					if ((paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id + 1 + x_count && Right_up == 1) || (paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id - 1 + x_count && Left_up == 1) || (paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id - 1 - x_count && Left_down == 1) || (paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id + 1 - x_count && Right_down == 1))
						break;
					if ((paths[i][j][paths[i][j].size() - 1].id % x_count == 0 && (paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id + 1 || paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id + 1 + x_count || paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id + 1 - x_count)) || (paths[i][j][paths[i][j].size() - 1].id % x_count == 1 && (paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id - 1 || paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id - 1 + x_count || paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id - 1 - x_count)))
						break;
					paths[i][j].insert(paths[i][j].end(), paths[i][m].begin(), paths[i][m].end());
					paths[i].erase(paths[i].begin() + m);
					m = 0; j = 0;
					break;
				}
			}
		}
		//合并不同斜率的路径
		for (int j = 0; j < paths[i].size(); j++) {
			int right_up = 0;
			int left_up = 0;
			int right_down = 0;
			int left_down = 0;
			for (int c = 0; c < paths[i].size(); c++) {
				for (int n = 1; n < paths[i][c].size(); n++) {
					if ((paths[i][j][0].id + 1 == paths[i][c][n - 1].id && paths[i][j][0].id + x_count == paths[i][c][n].id && paths[i][j][0].id % x_count != 0) ||
						(paths[i][j][0].id + 1 == paths[i][c][n].id && paths[i][j][0].id + x_count == paths[i][c][n - 1].id && paths[i][j][0].id % x_count != 0)) {
						right_up = 1;//如果与右上的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][0].id - 1 == paths[i][c][n - 1].id && paths[i][j][0].id + x_count == paths[i][c][n].id && paths[i][j][0].id % x_count != 1) ||
						(paths[i][j][0].id - 1 == paths[i][c][n].id && paths[i][j][0].id + x_count == paths[i][c][n - 1].id && paths[i][j][0].id % x_count != 1)) {
						left_up = 1;//如果与左上的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][0].id + 1 == paths[i][c][n - 1].id && paths[i][j][0].id - x_count == paths[i][c][n].id && paths[i][j][0].id % x_count != 0) ||
						(paths[i][j][0].id + 1 == paths[i][c][n].id && paths[i][j][0].id - x_count == paths[i][c][n - 1].id && paths[i][j][0].id % x_count != 0)) {
						right_down = 1;//如果与右下的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][0].id - 1 == paths[i][c][n - 1].id && paths[i][j][0].id - x_count == paths[i][c][n].id && paths[i][j][0].id % x_count != 1) ||
						(paths[i][j][0].id - 1 == paths[i][c][n].id && paths[i][j][0].id - x_count == paths[i][c][n - 1].id && paths[i][j][0].id % x_count != 1)) {
						left_down = 1;//如果与左下的端点连接，则会与已有的路径交叉
					}
				}
			}
			for (int m = 0; m < paths[i].size(); m++) {

				if (m == j) m++;
				if (m == paths[i].size()) break;
				//起点-起点
				if (abs(paths[i][m][1].id - paths[i][m][0].id) != abs(paths[i][j][0].id - paths[i][j][1].id) && paths[i][m][0].id == 2 * paths[i][j][0].id - paths[i][j][1].id) {
					if ((paths[i][m][0].id == paths[i][j][0].id + 1 + x_count && right_up == 1) || (paths[i][m][0].id == paths[i][j][0].id - 1 + x_count && left_up == 1) || (paths[i][m][0].id == paths[i][j][0].id - 1 - x_count && left_down == 1) || (paths[i][m][0].id == paths[i][j][0].id + 1 - x_count && right_down == 1))
						break;
					if ((paths[i][j][0].id % x_count == 0 && (paths[i][m][0].id == paths[i][j][0].id + 1 || paths[i][m][0].id == paths[i][j][0].id + 1 + x_count || paths[i][m][0].id == paths[i][j][0].id + 1 - x_count)) || (paths[i][j][0].id % x_count == 1 && (paths[i][m][0].id == paths[i][j][0].id - 1 || paths[i][m][0].id == paths[i][j][0].id - 1 + x_count || paths[i][m][0].id == paths[i][j][0].id - 1 - x_count)))
						break;
					reverse(paths[i][j].begin(), paths[i][j].end());//翻转
					paths[i][j].insert(paths[i][j].end(), paths[i][m].begin(), paths[i][m].end());
					paths[i].erase(paths[i].begin() + m);
					m = 0; j = 0;
					break;
				}
				//起点-终点
				if (abs(paths[i][m][paths[i][m].size() - 2].id - paths[i][m][paths[i][m].size() - 1].id) != abs(paths[i][j][0].id - paths[i][j][1].id) && paths[i][m][paths[i][m].size() - 1].id == 2 * paths[i][j][0].id - paths[i][j][1].id) {
					if ((paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id + 1 + x_count && right_up == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id - 1 + x_count && left_up == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id - 1 - x_count && left_down == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id + 1 - x_count && right_down == 1))
						break;
					if ((paths[i][j][0].id % x_count == 0 && (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id + 1 || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id + 1 + x_count || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id + 1 - x_count)) || (paths[i][j][0].id % x_count == 1 && (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id - 1 || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id - 1 + x_count || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id - 1 - x_count)))
						break;
					paths[i][m].insert(paths[i][m].end(), paths[i][j].begin(), paths[i][j].end());
					paths[i].erase(paths[i].begin() + j);
					m = 0; j = 0;
					break;
				}
			}
		}



		for (int j = 0; j < paths[i].size(); j++) {
			int right_up = 0;//判断该空点的斜连接是否可连接
			//可以连接为0，不可连接为1
			int left_up = 0;
			int right_down = 0;
			int left_down = 0;
			for (int c = 0; c < paths[i].size(); c++) {
				for (int n = 1; n < paths[i][c].size(); n++) {
					if ((paths[i][j][paths[i][j].size() - 1].id + 1 == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id + x_count == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 0) ||
						(paths[i][j][paths[i][j].size() - 1].id + 1 == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id + x_count == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 0)) {
						right_up = 1;//如果与右上的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][paths[i][j].size() - 1].id - 1 == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id + x_count == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 1) ||
						(paths[i][j][paths[i][j].size() - 1].id - 1 == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id + x_count == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 1)) {
						left_up = 1;//如果与左上的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][paths[i][j].size() - 1].id + 1 == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id - x_count == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 0) ||
						(paths[i][j][paths[i][j].size() - 1].id + 1 == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id - x_count == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 0)) {
						right_down = 1;//如果与右下的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][paths[i][j].size() - 1].id - 1 == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id - x_count == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 1) ||
						(paths[i][j][paths[i][j].size() - 1].id - 1 == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id - x_count == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 1)) {
						left_down = 1;//如果与左下的端点连接，则会与已有的路径交叉
					}
				}
			}
			for (int m = 0; m < paths[i].size(); m++) {



				if (m == j) m++;
				if (m == paths[i].size()) break;
				////终点-终点
				if (abs(paths[i][m][paths[i][m].size() - 2].id - paths[i][m][paths[i][m].size() - 1].id) != abs(paths[i][j][paths[i][j].size() - 1].id - paths[i][j][paths[i][j].size() - 2].id) && paths[i][m][paths[i][m].size() - 1].id == 2 * paths[i][j][paths[i][j].size() - 1].id - paths[i][j][paths[i][j].size() - 2].id) {
					if ((paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id + 1 + x_count && right_up == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id - 1 + x_count && left_up == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id - 1 - x_count && left_down == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id + 1 - x_count && right_down == 1))
						break;
					if ((paths[i][j][paths[i][j].size() - 1].id % x_count == 0 && (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id + 1 || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id + 1 + x_count || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id + 1 - x_count)) || (paths[i][j][paths[i][j].size() - 1].id % x_count == 1 && (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id - 1 || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id - 1 + x_count || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id - 1 - x_count)))
						break;
					reverse(paths[i][m].begin(), paths[i][m].end());//翻转
					paths[i][j].insert(paths[i][j].end(), paths[i][m].begin(), paths[i][m].end());
					paths[i].erase(paths[i].begin() + m);
					j = 0; m = 0;
					break;
				}
				////终点-起点
				if (abs(paths[i][m][1].id - paths[i][m][0].id) != abs(paths[i][j][paths[i][j].size() - 1].id - paths[i][j][paths[i][j].size() - 2].id) && paths[i][m][0].id == 2 * paths[i][j][paths[i][j].size() - 1].id - paths[i][j][paths[i][j].size() - 2].id) {
					if ((paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id + 1 + x_count && right_up == 1) || (paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id - 1 + x_count && left_up == 1) || (paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id - 1 - x_count && left_down == 1) || (paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id + 1 - x_count && right_down == 1))
						break;
					if ((paths[i][j][paths[i][j].size() - 1].id % x_count == 0 && (paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id + 1 || paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id + 1 + x_count || paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id + 1 - x_count)) || (paths[i][j][paths[i][j].size() - 1].id % x_count == 1 && (paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id - 1 || paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id - 1 + x_count || paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id - 1 - x_count)))
						break;
					paths[i][j].insert(paths[i][j].end(), paths[i][m].begin(), paths[i][m].end());
					paths[i].erase(paths[i].begin() + m);
					m = 0; j = 0;
					break;
				}
			}
		}
		

		//OptimizePoints();
		///////////////////////////////////////////////////
		for (int j = 0; j < model[i].size(); j++) {
		if (model[i][j].b == 0) {
			Right_up = 0;//判断该空点的斜连接是否可连接，可以连接为0，不可连接为1
			Left_up = 0;
			Right_down = 0;
			Left_down = 0;
			for (int m = 0; m < paths[i].size(); m++) {
				for (int n = 1; n < paths[i][m].size(); n++) {
					if ((model[i][j].id + 1 == paths[i][m][n - 1].id && model[i][j].id + x_count == paths[i][m][n].id && model[i][j].id % x_count != 0) ||
						(model[i][j].id + 1 == paths[i][m][n].id && model[i][j].id + x_count == paths[i][m][n - 1].id && model[i][j].id % x_count != 0)) {
						Right_up = 1;//如果与右上的端点连接，则会与已有的路径交叉
					}
					if ((model[i][j].id - 1 == paths[i][m][n - 1].id && model[i][j].id + x_count == paths[i][m][n].id && model[i][j].id % x_count != 1) ||
						(model[i][j].id - 1 == paths[i][m][n].id && model[i][j].id + x_count == paths[i][m][n - 1].id && model[i][j].id % x_count != 1)) {
						Left_up = 1;//如果与左上的端点连接，则会与已有的路径交叉
					}
					if ((model[i][j].id + 1 == paths[i][m][n - 1].id && model[i][j].id - x_count == paths[i][m][n].id && model[i][j].id % x_count != 0) ||
						(model[i][j].id + 1 == paths[i][m][n].id && model[i][j].id - x_count == paths[i][m][n - 1].id && model[i][j].id % x_count != 0)) {
						Right_down = 1;//如果与右下的端点连接，则会与已有的路径交叉
					}
					if ((model[i][j].id - 1 == paths[i][m][n - 1].id && model[i][j].id - x_count == paths[i][m][n].id && model[i][j].id % x_count != 1) ||
						(model[i][j].id - 1 == paths[i][m][n].id && model[i][j].id - x_count == paths[i][m][n - 1].id && model[i][j].id % x_count != 1)) {
						Left_down = 1;//如果与左下的端点连接，则会与已有的路径交叉
					}
				}
			}

			for (int m = 0; m < paths[i].size(); m++) {

				//右-终点
				if (paths[i][m][paths[i][m].size() - 1].id == model[i][j].id + 1 && model[i][j].id % x_count != 0) {
					model[i][j].b++;
					s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					paths[i][m].push_back(s);
					j = 0;
					break;//寻找到连接的路径，退出该空点的循环
				}//上-终点
				else if (paths[i][m][paths[i][m].size() - 1].id == model[i][j].id + x_count) {
					model[i][j].b++;
					s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					paths[i][m].push_back(s);
					j = 0;
					break;
				}//左-终点
				else if (paths[i][m][paths[i][m].size() - 1].id == model[i][j].id - 1 && model[i][j].id % x_count != 1) {
					model[i][j].b++;
					s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					paths[i][m].push_back(s);
					j = 0;
					break;
				}//下-终点
				else if (paths[i][m][paths[i][m].size() - 1].id == model[i][j].id - x_count) {
					model[i][j].b++;
					s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					paths[i][m].push_back(s);
					j = 0;
					break;
				}//右上-终点
				else if (paths[i][m][paths[i][m].size() - 1].id == model[i][j].id + 1 + x_count && model[i][j].id % x_count != 0 && Right_up == 0) {
					model[i][j].b++;
					s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					paths[i][m].push_back(s);
					j = 0;
					break;
				}//左上-终点
				else if (paths[i][m][paths[i][m].size() - 1].id == model[i][j].id - 1 + x_count && model[i][j].id % x_count != 1 && Left_up == 0) {
					model[i][j].b++;
					s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					paths[i][m].push_back(s);
					j = 0;
					break;
				}//左下-终点
				else if (paths[i][m][paths[i][m].size() - 1].id == model[i][j].id - 1 - x_count && model[i][j].id % x_count != 1 && Left_down == 0) {
					model[i][j].b++;
					s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					paths[i][m].push_back(s);
					j = 0;
					break;
				}//右下-终点
				else if (paths[i][m][paths[i][m].size() - 1].id == model[i][j].id + 1 - x_count && model[i][j].id % x_count != 0 && Right_down == 0) {
					model[i][j].b++;
					s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					paths[i][m].push_back(s);
					j = 0;
					break;
				}//右-起点
				else if (paths[i][m][0].id == model[i][j].id + 1 && model[i][j].id % x_count != 0) {
					model[i][j].b++;
					s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					paths[i][m].push_back(s);
					for (int f = paths[i][m].size() - 1; f > 0; f--) {
						s = paths[i][m][f];
						paths[i][m][f] = paths[i][m][f - 1];
						paths[i][m][f - 1] = s;
					}
					j = 0;
					break;
				}//上-起点
				else if (paths[i][m][0].id == model[i][j].id + x_count) {
					model[i][j].b++;
					s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					paths[i][m].push_back(s);
					for (int f = paths[i][m].size() - 1; f > 0; f--) {
						s = paths[i][m][f];
						paths[i][m][f] = paths[i][m][f - 1];
						paths[i][m][f - 1] = s;
					}
					j = 0;
					break;
				}//左-起点
				else if (paths[i][m][0].id == model[i][j].id - 1 && model[i][j].id % x_count != 1) {
					model[i][j].b++;
					s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					paths[i][m].push_back(s);
					for (int f = paths[i][m].size() - 1; f > 0; f--) {
						s = paths[i][m][f];
						paths[i][m][f] = paths[i][m][f - 1];
						paths[i][m][f - 1] = s;
					}
					j = 0;
					break;
				}//下-起点
				else if (paths[i][m][0].id == model[i][j].id - x_count) {
					model[i][j].b++;
					s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					paths[i][m].push_back(s);
					for (int f = paths[i][m].size() - 1; f > 0; f--) {
						s = paths[i][m][f];
						paths[i][m][f] = paths[i][m][f - 1];
						paths[i][m][f - 1] = s;
					}
					j = 0;
					break;
				}//右上-起点
				else if (paths[i][m][0].id == model[i][j].id + 1 + x_count && model[i][j].id % x_count != 0 && Right_up == 0) {
					model[i][j].b++;
					s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					paths[i][m].push_back(s);
					for (int f = paths[i][m].size() - 1; f > 0; f--) {
						s = paths[i][m][f];
						paths[i][m][f] = paths[i][m][f - 1];
						paths[i][m][f - 1] = s;
					}
					j = 0;
					break;
				}//左上-起点
				else if (paths[i][m][0].id == model[i][j].id - 1 + x_count && model[i][j].id % x_count != 1 && Left_up == 0) {
					model[i][j].b++;
					s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					paths[i][m].push_back(s);
					for (int f = paths[i][m].size() - 1; f > 0; f--) {
						s = paths[i][m][f];
						paths[i][m][f] = paths[i][m][f - 1];
						paths[i][m][f - 1] = s;
					}
					j = 0;
					break;
				}//左下-起点
				else if (paths[i][m][0].id == model[i][j].id - 1 - x_count && model[i][j].id % x_count != 1 && Left_down == 0) {
					model[i][j].b++;
					s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					paths[i][m].push_back(s);
					for (int f = paths[i][m].size() - 1; f > 0; f--) {
						s = paths[i][m][f];
						paths[i][m][f] = paths[i][m][f - 1];
						paths[i][m][f - 1] = s;
					}
					j = 0;
					break;
				}//右下-起点
				else if (paths[i][m][0].id == model[i][j].id + 1 - x_count && model[i][j].id % x_count != 0 && Right_down == 0) {
					model[i][j].b++;
					s.id = model[i][j].id; s.x = model[i][j].x; s.y = model[i][j].y; s.z = model[i][j].z; s.direction = model[i][j].direction;
					paths[i][m].push_back(s);
					for (int f = paths[i][m].size() - 1; f > 0; f--) {
						s = paths[i][m][f];
						paths[i][m][f] = paths[i][m][f - 1];
						paths[i][m][f - 1] = s;
					}
					j = 0;
					break;
				}


			}
		}
		}
		//最后合并
		for (int j = 0; j < paths[i].size(); j++) {
			Right_up = 0;
			Left_up = 0;
			Right_down = 0;
			Left_down = 0;
			for (int c = 0; c < paths[i].size(); c++) {
				for (int n = 1; n < paths[i][c].size(); n++) {
					if ((paths[i][j][0].id + 1 == paths[i][c][n - 1].id && paths[i][j][0].id + x_count == paths[i][c][n].id && paths[i][j][0].id % x_count != 0) ||
						(paths[i][j][0].id + 1 == paths[i][c][n].id && paths[i][j][0].id + x_count == paths[i][c][n - 1].id && paths[i][j][0].id % x_count != 0)) {
						Right_up = 1;//如果与右上的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][0].id - 1 == paths[i][c][n - 1].id && paths[i][j][0].id + x_count == paths[i][c][n].id && paths[i][j][0].id % x_count != 1) ||
						(paths[i][j][0].id - 1 == paths[i][c][n].id && paths[i][j][0].id + x_count == paths[i][c][n - 1].id && paths[i][j][0].id % x_count != 1)) {
						Left_up = 1;//如果与左上的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][0].id + 1 == paths[i][c][n - 1].id && paths[i][j][0].id - x_count == paths[i][c][n].id && paths[i][j][0].id % x_count != 0) ||
						(paths[i][j][0].id + 1 == paths[i][c][n].id && paths[i][j][0].id - x_count == paths[i][c][n - 1].id && paths[i][j][0].id % x_count != 0)) {
						Right_down = 1;//如果与右下的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][0].id - 1 == paths[i][c][n - 1].id && paths[i][j][0].id - x_count == paths[i][c][n].id && paths[i][j][0].id % x_count != 1) ||
						(paths[i][j][0].id - 1 == paths[i][c][n].id && paths[i][j][0].id - x_count == paths[i][c][n - 1].id && paths[i][j][0].id % x_count != 1)) {
						Left_down = 1;//如果与左下的端点连接，则会与已有的路径交叉
					}
				}
			}
			for (int m = 0; m < paths[i].size(); m++) {

				if (m == j) m++;
				if (m == paths[i].size()) break;
				//起点-起点
				if ((abs(paths[i][m][0].id - paths[i][j][0].id) == 1|| abs(paths[i][m][0].id - paths[i][j][0].id) == x_count|| abs(paths[i][m][0].id - paths[i][j][0].id) == x_count-1|| abs(paths[i][m][0].id - paths[i][j][0].id) == x_count+1)) {
					if ((paths[i][m][0].id == paths[i][j][0].id + 1 + x_count && Right_up == 1) || (paths[i][m][0].id == paths[i][j][0].id - 1 + x_count && Left_up == 1) || (paths[i][m][0].id == paths[i][j][0].id - 1 - x_count && Left_down == 1) || (paths[i][m][0].id == paths[i][j][0].id + 1 - x_count && Right_down == 1))
						break;
					if ((paths[i][j][0].id % x_count == 0&&(paths[i][m][0].id== paths[i][j][0].id+1|| paths[i][m][0].id == paths[i][j][0].id + 1+x_count|| paths[i][m][0].id == paths[i][j][0].id + 1-x_count) )||( paths[i][j][0].id % x_count == 1&&(paths[i][m][0].id == paths[i][j][0].id - 1|| paths[i][m][0].id == paths[i][j][0].id - 1+x_count|| paths[i][m][0].id == paths[i][j][0].id - 1-x_count)))
						break;
					reverse(paths[i][j].begin(), paths[i][j].end());//翻转
					paths[i][j].insert(paths[i][j].end(), paths[i][m].begin(), paths[i][m].end());
					paths[i].erase(paths[i].begin() + m);
					m = 0; j = 0;
					break;
				}
				//起点-终点
				if ( (abs(paths[i][m][paths[i][m].size() - 1].id - paths[i][j][0].id) == 1 || abs(paths[i][m][paths[i][m].size() - 1].id - paths[i][j][0].id) == x_count || abs(paths[i][m][paths[i][m].size() - 1].id - paths[i][j][0].id) == x_count - 1 || abs(paths[i][m][paths[i][m].size() - 1].id - paths[i][j][0].id) == x_count + 1)) {
					if ((paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id + 1 + x_count && Right_up == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id - 1 + x_count && Left_up == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id - 1 - x_count && Left_down == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id + 1 - x_count && Right_down == 1))
						break;
					if ((paths[i][j][0].id % x_count == 0 && (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id + 1 || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id + 1 + x_count || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id + 1 - x_count)) || (paths[i][j][0].id % x_count == 1 && (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id - 1 || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id - 1 + x_count || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][0].id - 1 - x_count)))
						break;
					paths[i][m].insert(paths[i][m].end(), paths[i][j].begin(), paths[i][j].end());
					paths[i].erase(paths[i].begin() + j);
					m = 0; j = 0;
					break;
				}
			}
		}



		for (int j = 0; j < paths[i].size(); j++) {
			Right_up = 0;//判断该空点的斜连接是否可连接
			//可以连接为0，不可连接为1
			Left_up = 0;
			Right_down = 0;
			Left_down = 0;
			for (int c = 0; c < paths[i].size(); c++) {
				for (int n = 1; n < paths[i][c].size(); n++) {
					if ((paths[i][j][paths[i][j].size() - 1].id + 1 == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id + x_count == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 0) ||
						(paths[i][j][paths[i][j].size() - 1].id + 1 == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id + x_count == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 0)) {
						Right_up = 1;//如果与右上的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][paths[i][j].size() - 1].id - 1 == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id + x_count == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 1) ||
						(paths[i][j][paths[i][j].size() - 1].id - 1 == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id + x_count == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 1)) {
						Left_up = 1;//如果与左上的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][paths[i][j].size() - 1].id + 1 == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id - x_count == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 0) ||
						(paths[i][j][paths[i][j].size() - 1].id + 1 == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id - x_count == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 0)) {
						Right_down = 1;//如果与右下的端点连接，则会与已有的路径交叉
					}
					if ((paths[i][j][paths[i][j].size() - 1].id - 1 == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id - x_count == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 1) ||
						(paths[i][j][paths[i][j].size() - 1].id - 1 == paths[i][c][n].id && paths[i][j][paths[i][j].size() - 1].id - x_count == paths[i][c][n - 1].id && paths[i][j][paths[i][j].size() - 1].id % x_count != 1)) {
						Left_down = 1;//如果与左下的端点连接，则会与已有的路径交叉
					}
				}
			}
			for (int m = 0; m < paths[i].size(); m++) {



				if (m == j) m++;
				if (m == paths[i].size()) break;
				////终点-终点
				if ((abs(paths[i][m][paths[i][m].size() - 1].id - paths[i][j][paths[i][j].size() - 1].id) == 1 || abs(paths[i][m][paths[i][m].size() - 1].id - paths[i][j][paths[i][j].size() - 1].id) == x_count || abs(paths[i][m][paths[i][m].size() - 1].id - paths[i][j][paths[i][j].size() - 1].id) == x_count - 1 || abs(paths[i][m][paths[i][m].size() - 1].id - paths[i][j][paths[i][j].size() - 1].id) == x_count + 1)) {
					if ((paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id + 1 + x_count && Right_up == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id - 1 + x_count && Left_up == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id - 1 - x_count && Left_down == 1) || (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id + 1 - x_count && Right_down == 1))
						break;
					if ((paths[i][j][paths[i][j].size() - 1].id % x_count == 0 && (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id + 1 || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id + 1 + x_count || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id + 1 - x_count)) || (paths[i][j][paths[i][j].size() - 1].id % x_count == 1 && (paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id - 1 || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id - 1 + x_count || paths[i][m][paths[i][m].size() - 1].id == paths[i][j][paths[i][j].size() - 1].id - 1 - x_count)))
						break;
					reverse(paths[i][m].begin(), paths[i][m].end());//翻转
					paths[i][j].insert(paths[i][j].end(), paths[i][m].begin(), paths[i][m].end());
					paths[i].erase(paths[i].begin() + m);
					j = 0; m = 0;
					break;
				}
				////终点-起点
				if ((abs(paths[i][m][0].id - paths[i][j][paths[i][j].size() - 1].id) == 1 || abs(paths[i][m][0].id - paths[i][j][paths[i][j].size() - 1].id) == x_count || abs(paths[i][m][0].id - paths[i][j][paths[i][j].size() - 1].id) == x_count - 1 || abs(paths[i][m][0].id - paths[i][j][paths[i][j].size() - 1].id) == x_count + 1)) {
					if ((paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id + 1 + x_count && Right_up == 1) || (paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id - 1 + x_count && Left_up == 1) || (paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id - 1 - x_count && Left_down == 1) || (paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id + 1 - x_count && Right_down == 1))
						break;
					if ((paths[i][j][paths[i][j].size() - 1].id % x_count == 0 && (paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id + 1 || paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id + 1 + x_count || paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id + 1 - x_count)) || (paths[i][j][paths[i][j].size() - 1].id % x_count == 1 && (paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id - 1 || paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id - 1 + x_count || paths[i][m][0].id == paths[i][j][paths[i][j].size() - 1].id - 1 - x_count)))
						break;
					paths[i][j].insert(paths[i][j].end(), paths[i][m].begin(), paths[i][m].end());
					paths[i].erase(paths[i].begin() + m);
					m = 0; j = 0;
					break;
				}
			}
		}

		printf("%.2lf%%\r", i * 100.0 / paths.size());
	}
	
}

void StressOptimize() {
	//根据应力对应力小的区域进行过滤
	double v = 0.3;//过滤百分比,gong = 0.05;ls = ? ; qiao = 0.1
	double a;
	bool flag;
	element s;
	vector<element>q;//将每个离散单元存入q中，便于计算改阈值下的应力大小
	//for (int i = 0; i < point.size(); i++) {//将每个单元的数据存储到q中
	for (int j = 0; j < point[0].size(); j++) {
		for (int m = 0; m < point[0][j].size(); m++) {
			if (point[0][j][m].b == 1)
				q.push_back(point[0][j][m]);
		}
	}
	//}
	printf("对缓存进行排序：\n");
	for (int j = 0; j < q.size() - 1; j++) {//对q进行从小到达排序
		flag = false;
		for (int m = q.size() - 1; m > j; m--) {
			if (abs(q[m - 1].s_max) > abs(q[m].s_max)) {
				s = q[m - 1];
				q[m - 1] = q[m];
				q[m] = s;
				flag = true;
			}
		}
		if (flag == false) break;
		printf("%.2lf%%\r", j * 100.0 / (q.size() - 1));
	}
	a = abs(q[0].s_max) + (abs(q[q.size() - 1].s_max) - abs(q[0].s_max)) * v;
	printf("优化应力小区域：\n");
	//方法一
	//方法一（迭代，使用广度优先遍历算法）
	//for (int i = 0; i < model.size(); i++) {
	//	for (int j = 0; j < model[i].size(); j++) {
	//		if (abs(model[i][j].s_max) <= a) {
	//			model[i][j].direction = -1;
	//			que.push(model[i][j]);//存入队列
	//		}
	//	}
	//	while (!que.empty()) {
	//		points front = que.front();//第一个出队
	//		que.pop();//删除第一个元素
	//		double count = 0, sum = 0;
	//		//遍历当前空单元的八个邻接单元
	//		for (int m = 0; m < model[i].size(); m++) {
	//			if (model[i][m].direction != -1) {
	//				if (model[i][m].id == front.id + x_count || model[i][m].id == front.id - x_count ||
	//					(model[i][m].id == front.id + x_count + 1 && front.id % x_count != 0) || (model[i][m].id == front.id + 1 && front.id % x_count != 0) ||
	//					(model[i][m].id == front.id - x_count + 1 && front.id % x_count != 0) || (model[i][m].id == front.id - x_count - 1 && front.id % x_count != 1) ||
	//					(model[i][m].id == front.id - 1 && front.id % x_count != 1) || (model[i][m].id == front.id - 1 + x_count && front.id % x_count != 1)) {
	//					sum += model[i][m].direction;
	//					count++;
	//				}
	//			}
	//		}
	//		if (count == 0) {
	//			que.push(front);
	//		}
	//		else {
	//			for (int m = 0; m < model[i].size(); m++) {
	//				if (model[i][m].id == front.id)
	//					model[i][m].direction = sum / count;
	//			}
	//		}
	//	}
	//	printf("%.2lf%%\r", i * 100.0 / model.size());
	//}
	//方法二
	for (int i = 0; i < model.size(); i++) {
		for (int j = 0; j < model[i].size(); j++) {
			if (abs(model[i][j].s_max) <= a) {
				model[i][j].direction = -1;
			}
		}
		for (int j = 0; j < model[i].size(); j++) {
			if (model[i][j].direction != -1) {
				for (int m = 0; m < model[i].size(); m++) {
					if (model[i][m].direction == -1) {
						if (model[i][m].id == model[i][j].id + 1 || model[i][m].id == model[i][j].id - 1) {
							model[i][m].direction = model[i][j].direction;
					
						}
						else if ((model[i][m].id == model[i][j].id + 1 - x_count || model[i][m].id == model[i][j].id + 1 || model[i][m].id == model[i][j].id + 1 + x_count) && model[i][j].id % x_count != 0) {
							model[i][m].direction = model[i][j].direction;
							
						}
						else if ((model[i][m].id == model[i][j].id - 1 - x_count || model[i][m].id == model[i][j].id - 1 || model[i][m].id == model[i][j].id - 1 + x_count) && model[i][j].id % x_count != 1) {
							model[i][m].direction = model[i][j].direction;
							
						}

					}
				}
			}
			if (j == model[i].size()-1) {
				for (int n = 0; n < model[i].size(); n++) {
					if (model[i][n].direction == -1) {
						j = 0;
						break;
					}
				}
			}
			
		}
		printf("%.2lf%%\r", i * 100.0 / model.size());
	}


}

void StressPath() {
	points s;//工具
	path n;
	bool flag;
	
	//将所有bool为1的单元存入model中
	for (int i = 0; i < point.size(); i++) {
		for (int j = 0; j < point[0].size(); j++) {
			for (int m = 0; m < point[0][0].size(); m++) {
				if (point[i][j][m].b == 1) {
					s.id = point[i][j][m].id;
					s.x = point[i][j][m].x;
					s.y = point[i][j][m].y;
					s.z = point[i][j][m].z;
					s.s_max = point[i][j][m].s_max;
					s.direction = point[i][j][m].direction;
					model_buffer.push_back(s);
				}

			}
		}
		model.push_back(model_buffer);
		model_buffer.clear();
	}
	//按照应力（abs）的大小对其进行排序
	//采用冒泡排序;
	printf("根据应力大小进行排序：\n");
	for (int i = 0; i < model.size(); i++) {
		
		for (int j = 0; j < model[i].size() - 1; j++) {
			flag = false;
			for (int m = model[i].size() - 1; m > j; m--) {
				if (abs(model[i][m - 1].s_max) < abs(model[i][m].s_max)) {
					s = model[i][m - 1];
					model[i][m - 1] = model[i][m];
					model[i][m] = s;
					flag = true;
				}
			}
			if (flag == false) break;
		}
		printf("%.2lf%%\r", i * 100.0 / model.size());
	}

	//StressOptimize();//对应力小的区域进行优化


	//生成应力路径
	for (int i = 0; i < model.size(); i++) {
		for (int j = 0; j < model[i].size(); j++) {
			if (model[i][j].b == 0) {
				//model[i][j].b++;//增加该单元的度数
				s = model[i][j];
				n.id = s.id;
				n.x = s.x;
				n.y = s.y;
				n.z = s.z;
				n.direction = s.direction;
				paths_buffer1.push_back(n);
				VectorInsert(i, j);
				paths_buffer2.push_back(paths_buffer1);
				paths_buffer1.clear();
			}
			
		}
		printf("%.2lf%%\r", i * 100.0 / model.size());
		paths.push_back(paths_buffer2);
		paths_buffer2.clear();
	}
	//删除只有一个单元的路径
	for (int i = 0; i < paths.size(); i++) {
		for (int j = 0; j < paths[i].size(); j++) {
			if (paths[i][j].size() == 1) {
				paths[i].erase(paths[i].begin() + j);
				j--;//每删除一个元素，其后面的地址会向前进一个单位，所以j需要-1
			}
		}
	}
	//优化路径空点
	printf("优化空点......\n");
	OptimizePoints();
	//连续路径优化
	printf("连续路径优化......\n");

	PathConnect();
	//优化G-code路径
	printf("连续G-code路径......\n");
	/*for (int i = 0; i < paths.size(); i++) {
		for (int j = 0; j < paths[i].size(); j++) {
			if (paths[i][j].size() > 2) {
				for (int m = 2; m < paths[i][j].size(); m++) {
					if (((paths[i][j][m].x==paths[i][j][m-1].x&&paths[i][j][m-1].x==paths[i][j][m-2].x)||(paths[i][j][m].y == paths[i][j][m - 1].y && paths[i][j][m - 1].y == paths[i][j][m - 2].y))
					) {
						paths[i][j].erase(paths[i][j].begin() + m - 1);
						m = m - 1;
					}
				}
			}

		}
	}*/

	

}

void GcodePrint() {
	FILE* fp;
	errno_t err;     //判断此文件流是否存在 存在返回1
	err = fopen_s(&fp, "ls2.gcode", "a"); //若return 1 , 则将指向这个文件的文件流给
	double t1 = 0.03326;//层厚0.2，丝宽0.4,挤出轮步进
	double t2 = 0.04756;//层厚0.2，丝宽0.56
	double t0 = t1 + (t2 - t1) / 2;
	double E = 0;
	double r;//回抽
	int L=-50;//偏移量
	fprintf(fp, ";FLAVOR:Marlin\n");
	fprintf(fp, ";Generated with Cura_SteamEngine 4.10.0\n");
	fprintf(fp, "M140 S50\n");
	fprintf(fp, "M105\n");
	fprintf(fp, "M190 S50\n");
	fprintf(fp, "M104 S210\n");
	fprintf(fp, "M105\n");
	fprintf(fp, "M109 S210\n");
	fprintf(fp, "M82 ;absolute extrusion mode\n");
	fprintf(fp, "M201 X500.00 Y500.00 Z100.00 E5000.00 ;Setup machine max acceleration\n");
	fprintf(fp, "M203 X500.00 Y500.00 Z10.00 E50.00 ;Setup machine max feedrate\n");
	fprintf(fp, "M204 P500.00 R1000.00 T500.00 ;Setup Print/Retract/Travel acceleration\n");
	fprintf(fp, "M205 X8.00 Y8.00 Z0.40 E5.00 ;Setup Jerk\n");
	fprintf(fp, "M220 S100 ;Reset Feedrate\n");
	fprintf(fp, "M221 S100 ;Reset Flowrate\n");

	fprintf(fp, "G28 ;Home\n");

	fprintf(fp, "G92 E0\n");
	fprintf(fp, "G92 E0\n");
	fprintf(fp, "G1 F2700 E-5\n");
	fprintf(fp, "M107\n");
	

	
	
	//模型第一层
	
	//fprintf(fp, "G0 F9000 X%f Y%f Z%f\n", paths[0][0][0].x+L, paths[0][0][0].y + L, paths[0][0][0].z + 0.1);//模型第一层与底座间隔0.22
	//fprintf(fp, "G1 F4800 E0.00000\n");
	//
	//fprintf(fp, ";TYPE:WALL-OUTPE\n");


	//for (int j = 0; j < path_w.size(); j++) {
	//	r = E - 2;
	//	fprintf(fp, "G1 F4800 E%f\n", r);//回抽一定距离，避免拉丝
	//	fprintf(fp, "G0 X%f Y%f\n", path_w[j][0].x + L, path_w[j][0].y + L);
	//	fprintf(fp, "G1 F4800 E%f\n", E);
	//	for (int m = 1; m < path_w[j].size(); m++) {
	//		fprintf(fp, "G1 F1200 X%f Y%f E%f\n", path_w[j][m].x + L, path_w[j][m].y + L, E += distance(path_w[j][m - 1], path_w[j][m]) * t2);
	//	}
	//}


	//fprintf(fp, ";TYPE:FILL\n");
	//for (int j = 0; j < paths[0].size(); j++) {
	//	r = E - 2;
	//	fprintf(fp, "G1 F4800 E%f\n",r);//回抽一定距离，避免拉丝
	//	fprintf(fp, "G0 X%f Y%f\n", paths[0][j][0].x + L, paths[0][j][0].y + L);
	//	fprintf(fp, "G1 F4800 E%f\n", E);
	//	for (int m = 1; m < paths[0][j].size(); m++) {
	//		if (paths[0][j][m - 1].x != paths[0][j][m].x && paths[0][j][m - 1].y != paths[0][j][m].y) {
	//			fprintf(fp, "G1 F1200 X%f Y%f E%f\n", paths[0][j][m].x + L, paths[0][j][m].y + L, E += distance(paths[0][j][m - 1], paths[0][j][m]) * t1);
	//		}
	//		else {
	//			fprintf(fp, "G1 F1200 X%f Y%f E%f\n", paths[0][j][m].x + L, paths[0][j][m].y + L, E += distance(paths[0][j][m - 1], paths[0][j][m]) * t2);
	//		}
	//	}
	//}

	for (int i = 5; i < paths.size(); i++) {//从第二层开始
		//fprintf(fp, "G0 F9000 X%f Y%f Z%f\n", path_w[i][0][0].x + L, path_w[i][0][0].y + L, paths[i][0][0].z + 0.1);
		//fprintf(fp, ";TYPE:WALL-OUTPE\n");
		//for (int j = 0; j < path_w[i].size(); j++) {

		//	r = E - 2;
		//	fprintf(fp, "G1 F4800 E%f\n", r);//回抽一定距离，避免拉丝
		//	fprintf(fp, "G0 X%f Y%f\n", path_w[i][j][0].x + L, path_w[i][j][0].y + L);
		//	fprintf(fp, "G1 F4800 E%f\n", E);
		//	for (int m = 1; m < path_w[i][j].size(); m++) {
		//		fprintf(fp, "G1 F1200 X%f Y%f E%f\n", path_w[i][j][m].x + L, path_w[i][j][m].y + L, E += distance(path_w[i][j][m - 1], path_w[i][j][m]) * t2);
		//	}
		//}

		fprintf(fp, "G0 F9000 X%f Y%f Z%f\n", paths[i][0][0].x + L, paths[i][0][0].y + L, paths[i][0][0].z + 0.1);
		fprintf(fp, ";TYPE:FILL\n");
		for (int j = 0; j < paths[i].size(); j++) {

			r = E - 2;
			fprintf(fp, "G1 F4800 E%f\n", r);//回抽一定距离，避免拉丝
			fprintf(fp, "G0 X%f Y%f\n", paths[i][j][0].x + L, paths[i][j][0].y + L);
			fprintf(fp, "G1 F4800 E%f\n", E);

			for (int m = 1; m < paths[i][j].size(); m++) {
				if (paths[i][j][m - 1].x != paths[i][j][m].x && paths[i][j][m - 1].y != paths[i][j][m].y) {
					
					fprintf(fp, "G1 F3000 X%f Y%f E%f\n", paths[i][j][m].x + L, paths[i][j][m].y + L, E += distance(paths[i][j][m - 1], paths[i][j][m]) * t1);
				}
				else {
					
					fprintf(fp, "G1 F3000 X%f Y%f E%f\n", paths[i][j][m].x + L, paths[i][j][m].y + L, E += distance(paths[i][j][m - 1], paths[i][j][m]) * t2);
				}
			}
		}
	}
	fprintf(fp, "M140 S0\n");
	fprintf(fp, "M107\n");
	fprintf(fp, "G91\n");
	fprintf(fp, "G1 E-2 F2700\n");
	fprintf(fp, "G1 E-2 Z0.2 F2400 ;Retract and raise Z\n");
	fprintf(fp, "G1 X5 Y5 F3000 ;Wipe out\n");
	fprintf(fp, "G1 Z10 ;Raise Z more\n");
	fprintf(fp, "G90 ;Absolute positioning\n");

	fprintf(fp, "G1 X0 Y300 ;Present print\n");
	fprintf(fp, "M106 S0 ;Turn-off fan\n");
	fprintf(fp, "M104 S0 ;Turn-off hotend\n");
	fprintf(fp, "M140 S0 ;Turn-off bed\n");

	fprintf(fp, "M84 X Y E ;Disable all steppers but Z\n");

	fprintf(fp, "M82 ;absolute extrusion mode\n");
	fprintf(fp, "M104 S0\n");
	fclose(fp);
}

void AbaqusData() {
	FILE* a; 
	errno_t err;//判断此文件流是否存在 存在返回1
	err = fopen_s(&a, "Abaqus_a.txt", "a"); //若return 1 , 则将指向这个文件的文件流给
	FILE* b;    //判断此文件流是否存在 存在返回1
	err = fopen_s(&b, "Abaqus_b.txt", "a"); //若return 1 , 则将指向这个文件的文件流给
	FILE* c;
	err = fopen_s(&c, "Abaqus_c.txt", "a");
	float b1, b2, c1;
	for (int i = 0; i < 2; i++) {

		for (int j = 0; j < point[i].size(); j++) {
			for (int m = 0; m < point[i][j].size(); m++) {
				fprintf(a, "%d\n", point[i][j][m].b);
				b1 = (double)180 * (pi - point[i][j][m].direction) / pi;
				if ((b1 >= 0 && b1 < 22.5) || (b1 >= 157.5 && b1 < 180)) b2 = 0;
				else if (b1 >= 22.5 && b1 < 67.5) b2 = 45;
				else if (b1 >= 67.5 && b1 < 112.5) b2 = 90;
				else if (b1 >= 112.5 && b1 < 157.5) b2 = 135;
				fprintf(b, "%.1f\n", b2);
				if (b1 >= 0 && b1 < 60) c1 = 30;
				else if (b1 >= 60 && b1 < 120) c1 = 90;
				else if (b1 >= 120 && b1 < 180) c1 = 150;
				fprintf(c, "%.1f\n", c1);
			}
		}
	}
	fclose(a);
	fclose(b);
	fclose(c);
}

void MatlabData() {
	int d = 0;
	int num;
	FILE* a1;
	errno_t err;//判断此文件流是否存在 存在返回1
	err = fopen_s(&a1, "Matlab_x.txt", "a"); //若return 1 , 则将指向这个文件的文件流给
	FILE* b1;    //判断此文件流是否存在 存在返回1
	err = fopen_s(&b1, "Matlab_y.txt", "a"); //若return 1 , 则将指向这个文件的文件流给
	num = paths[paths.size()-1][0].size();
	for (int i = 0; i < paths[paths.size() - 1].size(); i++) {
		if (paths[paths.size() - 1][i].size() > num)
			num = paths[paths.size() - 1][i].size();
	}

	for (int i = 0; i < paths[paths.size() - 1].size(); i++) {
		for (int j = 0; j < num; j++) {
			if (j < paths[paths.size() - 1][i].size()) {
				fprintf(a1, "%f\t", paths[paths.size() - 1][i][j].x);
				fprintf(b1, "%f\t", paths[paths.size() - 1][i][j].y);
			}
			else {
				fprintf(a1, "%f\t",0);
				fprintf(b1, "%f\t",0);
			}
			

		}
		fprintf(a1, "\n");
		fprintf(b1, "\n");
	}
	fclose(a1);
	fclose(b1);
}

void fangxiangData() {
	int d = 0;
	int num;
	FILE* a1;
	errno_t err;//判断此文件流是否存在 存在返回1
	err = fopen_s(&a1, "fx.txt", "a"); //若return 1 , 则将指向这个文件的文件流给
	
	for (int j = 0; j < model[0].size(); j++) {
		fprintf(a1, "%f %f %f %f\n", model[0][j].x, model[0][j].y,model[0][j].s_max, model[0][j].direction);
	}
	fclose(a1);

}


void double_path() {
	path p;
	int t = 0;//用于判断起点和终点
	paths_buffer1.clear();
	paths_buffer2.clear();
	paths_buffer3.clear();
	for (int i = 0; i < paths.size(); i++) {
		for (int j = 0; j < paths[i].size(); j++) {
			for (int n = 0; n < 2; n++) {
				if (n == 1) {
					reverse(paths[i][j].begin(), paths[i][j].end());
				}
				for (int m = 0; m < paths[i][j].size(); m++) {
					if (m == 0) {//起点
						if (paths[i][j][m + 1].id == paths[i][j][m].id - x_count) {//上->下
							p.x = paths[i][j][m].x - d1 * 0.5;
							p.y = paths[i][j][m].y + d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							p.x = paths[i][j][m].x - d1 * 0.5;
							p.y = paths[i][j][m].y;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if (paths[i][j][m + 1].id == paths[i][j][m].id + x_count) {//下->上
							p.x = paths[i][j][m].x + d1 * 0.5;
							p.y = paths[i][j][m].y - d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							p.x = paths[i][j][m].x + d1 * 0.5;
							p.y = paths[i][j][m].y;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if (paths[i][j][m + 1].id == paths[i][j][m].id + 1) {//左->右
							p.x = paths[i][j][m].x - d1 * 0.5;
							p.y = paths[i][j][m].y - d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							p.x = paths[i][j][m].x;
							p.y = paths[i][j][m].y - d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						
						}
						else if (paths[i][j][m + 1].id == paths[i][j][m].id - 1) {//右->左
							p.x = paths[i][j][m].x + d1 * 0.5;
							p.y = paths[i][j][m].y + d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							p.x = paths[i][j][m].x;
							p.y = paths[i][j][m].y + d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if (paths[i][j][m + 1].id == paths[i][j][m].id - x_count + 1) {//左上->右下
							p.x = paths[i][j][m].x - sqrt(2) * d2 * 0.5;
							p.y = paths[i][j][m].y;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							p.x = paths[i][j][m].x - (d2 * 0.5) * (sqrt(2) * 0.5);
							p.y = paths[i][j][m].y - (d2 * 0.5) * (sqrt(2) * 0.5);
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if (paths[i][j][m + 1].id == paths[i][j][m].id + x_count - 1) {//右下->左上
							p.x = paths[i][j][m].x + sqrt(2) * d2 * 0.5;
							p.y = paths[i][j][m].y;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							p.x = paths[i][j][m].x + (d2 * 0.5) * (sqrt(2) * 0.5);
							p.y = paths[i][j][m].y + (d2 * 0.5) * (sqrt(2) * 0.5);
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if (paths[i][j][m + 1].id == paths[i][j][m].id - x_count - 1) {//右上->左下
							p.x = paths[i][j][m].x;
							p.y = paths[i][j][m].y + sqrt(2) * d2 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							p.x = paths[i][j][m].x - (d2 * 0.5) * (sqrt(2) * 0.5);
							p.y = paths[i][j][m].y + (d2 * 0.5) * (sqrt(2) * 0.5);
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if (paths[i][j][m + 1].id == paths[i][j][m].id + x_count + 1) {//左下->右上
							p.x = paths[i][j][m].x;
							p.y = paths[i][j][m].y - sqrt(2) * d2 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							p.x = paths[i][j][m].x + (d2 * 0.5) * (sqrt(2) * 0.5);
							p.y = paths[i][j][m].y - (d2 * 0.5) * (sqrt(2) * 0.5);
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
					}
					else if (m == paths[i][j].size() - 1) {//终点
						if (paths[i][j][m - 1].id == paths[i][j][m].id + x_count) {//上-下
							p.x = paths[i][j][m].x - d1 * 0.5;
							p.y = paths[i][j][m].y;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							p.x = paths[i][j][m].x - d1 * 0.5;
							p.y = paths[i][j][m].y - d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						
						}
						else if (paths[i][j][m - 1].id == paths[i][j][m].id - x_count) {//下-上
							p.x = paths[i][j][m].x + d1 * 0.5;
							p.y = paths[i][j][m].y;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							p.x = paths[i][j][m].x + d1 * 0.5;
							p.y = paths[i][j][m].y + d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							
						}
						else if (paths[i][j][m - 1].id == paths[i][j][m].id - 1) {//左-右
							p.x = paths[i][j][m].x;
							p.y = paths[i][j][m].y - d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							p.x = paths[i][j][m].x + d1 * 0.5;
							p.y = paths[i][j][m].y - d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						
						}
						else if (paths[i][j][m - 1].id == paths[i][j][m].id + 1) {//右-左
							p.x = paths[i][j][m].x;
							p.y = paths[i][j][m].y + d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							p.x = paths[i][j][m].x - d1 * 0.5;
							p.y = paths[i][j][m].y + d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							
						}
						else if (paths[i][j][m - 1].id == paths[i][j][m].id + x_count - 1) {//左上-右下
							p.x = paths[i][j][m].x - (d2 * 0.5) * (sqrt(2) * 0.5);
							p.y = paths[i][j][m].y - (d2 * 0.5) * (sqrt(2) * 0.5);
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							p.x = paths[i][j][m].x;
							p.y = paths[i][j][m].y - sqrt(2) * d2 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							
						}
						else if (paths[i][j][m - 1].id == paths[i][j][m].id - x_count + 1) {//右下-左上
							p.x = paths[i][j][m].x + (d2 * 0.5) * (sqrt(2) * 0.5);
							p.y = paths[i][j][m].y + (d2 * 0.5) * (sqrt(2) * 0.5);
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							p.x = paths[i][j][m].x;
							p.y = paths[i][j][m].y + sqrt(2) * d2 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						
						}
						else if (paths[i][j][m - 1].id == paths[i][j][m].id + x_count + 1) {//右上-左下
							p.x = paths[i][j][m].x - (d2 * 0.5) * (sqrt(2) * 0.5);
							p.y = paths[i][j][m].y + (d2 * 0.5) * (sqrt(2) * 0.5);
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							p.x = paths[i][j][m].x - sqrt(2) * d2 * 0.5;
							p.y = paths[i][j][m].y;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							
						}
						else if (paths[i][j][m - 1].id == paths[i][j][m].id - x_count - 1) {//左下-右上
							p.x = paths[i][j][m].x + (d2 * 0.5) * (sqrt(2) * 0.5);
							p.y = paths[i][j][m].y - (d2 * 0.5) * (sqrt(2) * 0.5);
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
							p.x = paths[i][j][m].x + sqrt(2) * d2 * 0.5;
							p.y = paths[i][j][m].y;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						
						}
					}
					else {
						if ((paths[i][j][m - 1].id == paths[i][j][m].id + x_count && paths[i][j][m].id == paths[i][j][m + 1].id + x_count) || 
							(paths[i][j][m - 1].id == paths[i][j][m].id + x_count && paths[i][j][m].id == paths[i][j][m + 1].id + x_count - 1) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id + x_count && paths[i][j][m].id == paths[i][j][m + 1].id + x_count + 1) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id + x_count - 1 && paths[i][j][m].id == paths[i][j][m + 1].id + x_count) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id + x_count + 1 && paths[i][j][m].id == paths[i][j][m + 1].id + x_count)) { //case 1
							p.x = paths[i][j][m].x - d1 * 0.5;
							p.y = paths[i][j][m].y;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id - x_count && paths[i][j][m].id == paths[i][j][m + 1].id - x_count) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id - x_count + 1 && paths[i][j][m].id == paths[i][j][m + 1].id - x_count) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id - x_count - 1 && paths[i][j][m].id == paths[i][j][m + 1].id - x_count) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id - x_count && paths[i][j][m].id == paths[i][j][m + 1].id - x_count + 1) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id - x_count && paths[i][j][m].id == paths[i][j][m + 1].id - x_count - 1)) {//case 2
							p.x = paths[i][j][m].x + d1 * 0.5;
							p.y = paths[i][j][m].y;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id - 1 && paths[i][j][m].id == paths[i][j][m + 1].id - 1) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id - 1 && paths[i][j][m].id == paths[i][j][m + 1].id - 1 - x_count) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id - 1 && paths[i][j][m].id == paths[i][j][m + 1].id - 1 + x_count) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id - 1 + x_count && paths[i][j][m].id == paths[i][j][m + 1].id - 1) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id - 1 - x_count && paths[i][j][m].id == paths[i][j][m + 1].id - 1)) {//case 3
							p.x = paths[i][j][m].x;
							p.y = paths[i][j][m].y - d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id + 1 && paths[i][j][m].id == paths[i][j][m + 1].id + 1) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id + 1 + x_count && paths[i][j][m].id == paths[i][j][m + 1].id + 1) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id + 1 - x_count && paths[i][j][m].id == paths[i][j][m + 1].id + 1) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id + 1 && paths[i][j][m].id == paths[i][j][m + 1].id + 1 - x_count) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id + 1 && paths[i][j][m].id == paths[i][j][m + 1].id + 1 + x_count)) {//case 4
							p.x = paths[i][j][m].x;
							p.y = paths[i][j][m].y + d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if (paths[i][j][m - 1].id == paths[i][j][m].id + x_count - 1 && paths[i][j][m].id == paths[i][j][m + 1].id + x_count - 1) {//case 5
							p.x = paths[i][j][m].x - (d2 * 0.5) * (sqrt(2) * 0.5);
							p.y = paths[i][j][m].y - (d2 * 0.5) * (sqrt(2) * 0.5);
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if (paths[i][j][m - 1].id == paths[i][j][m].id - x_count + 1 && paths[i][j][m].id == paths[i][j][m + 1].id - x_count + 1) {//case 6
							p.x = paths[i][j][m].x + (d2 * 0.5) * (sqrt(2) * 0.5);
							p.y = paths[i][j][m].y + (d2 * 0.5) * (sqrt(2) * 0.5);
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if (paths[i][j][m - 1].id == paths[i][j][m].id - x_count - 1 && paths[i][j][m].id == paths[i][j][m + 1].id - x_count - 1) {//case 7
							p.x = paths[i][j][m].x + (d2 * 0.5) * (sqrt(2) * 0.5);
							p.y = paths[i][j][m].y - (d2 * 0.5) * (sqrt(2) * 0.5);
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if (paths[i][j][m - 1].id == paths[i][j][m].id + x_count + 1 && paths[i][j][m].id == paths[i][j][m + 1].id + x_count + 1) {//case 8
							p.x = paths[i][j][m].x - (d2 * 0.5) * (sqrt(2) * 0.5);
							p.y = paths[i][j][m].y + (d2 * 0.5) * (sqrt(2) * 0.5);
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id + x_count && paths[i][j][m].id == paths[i][j][m + 1].id - 1) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id - 1 && paths[i][j][m].id == paths[i][j][m + 1].id + x_count)) {//case 9
							p.x = paths[i][j][m].x - d1 * 0.5;
							p.y = paths[i][j][m].y - d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id + 1 && paths[i][j][m].id == paths[i][j][m + 1].id - x_count) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id - x_count && paths[i][j][m].id == paths[i][j][m + 1].id + 1)) {//case 10
							p.x = paths[i][j][m].x + d1 * 0.5;
							p.y = paths[i][j][m].y + d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id - 1 && paths[i][j][m].id == paths[i][j][m + 1].id - x_count) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id - x_count && paths[i][j][m].id == paths[i][j][m + 1].id - 1)) {//case 11
							p.x = paths[i][j][m].x + d1 * 0.5;
							p.y = paths[i][j][m].y - d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id + x_count && paths[i][j][m].id == paths[i][j][m + 1].id + 1) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id + 1 && paths[i][j][m].id == paths[i][j][m + 1].id + x_count)) {//case 12
							p.x = paths[i][j][m].x - d1 * 0.5;
							p.y = paths[i][j][m].y + d1 * 0.5;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id + x_count && paths[i][j][m].id == paths[i][j][m + 1].id - x_count - 1)||
							(paths[i][j][m - 1].id == paths[i][j][m].id - x_count - 1 && paths[i][j][m].id == paths[i][j][m + 1].id + x_count)) {//case 13
							p.x = paths[i][j][m].x - d1 * 0.5;
							p.y = paths[i][j][m].y - d1;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id + x_count +1 && paths[i][j][m].id == paths[i][j][m + 1].id - x_count) ||
							(paths[i][j][m - 1].id == paths[i][j][m].id - x_count && paths[i][j][m].id == paths[i][j][m + 1].id + x_count + 1)) {//case 14
							p.x = paths[i][j][m].x + d1 * 0.5;
							p.y = paths[i][j][m].y + d1;
							p.z = paths[i][j][m].z;
							paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id + x_count && paths[i][j][m].id == paths[i][j][m + 1].id - x_count + 1) ||
						(paths[i][j][m - 1].id == paths[i][j][m].id - x_count + 1 && paths[i][j][m].id == paths[i][j][m + 1].id + x_count)) {//case 15
						p.x = paths[i][j][m].x - d1 * 0.5;
						p.y = paths[i][j][m].y + d1;
						p.z = paths[i][j][m].z;
						paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id + x_count - 1 && paths[i][j][m].id == paths[i][j][m + 1].id - x_count) ||
						(paths[i][j][m - 1].id == paths[i][j][m].id - x_count && paths[i][j][m].id == paths[i][j][m + 1].id + x_count - 1)) {//case 16
						p.x = paths[i][j][m].x + d1 * 0.5;
						p.y = paths[i][j][m].y - d1;
						p.z = paths[i][j][m].z;
						paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id  - 1 && paths[i][j][m].id == paths[i][j][m + 1].id - x_count + 1) ||
						(paths[i][j][m - 1].id == paths[i][j][m].id - x_count + 1 && paths[i][j][m].id == paths[i][j][m + 1].id - 1)) {//case 17
						p.x = paths[i][j][m].x + d1;
						p.y = paths[i][j][m].y - d1 * 0.5;
						p.z = paths[i][j][m].z;
						paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id - 1 + x_count && paths[i][j][m].id == paths[i][j][m + 1].id  + 1) ||
						(paths[i][j][m - 1].id == paths[i][j][m].id + 1 && paths[i][j][m].id == paths[i][j][m + 1].id - 1 + x_count)) {//case 18
						p.x = paths[i][j][m].x - d1;
						p.y = paths[i][j][m].y + d1 * 0.5;
						p.z = paths[i][j][m].z;
						paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id + 1 + x_count && paths[i][j][m].id == paths[i][j][m + 1].id - 1) ||
						(paths[i][j][m - 1].id == paths[i][j][m].id - 1 && paths[i][j][m].id == paths[i][j][m + 1].id + 1 + x_count)) {//case 19
						p.x = paths[i][j][m].x - d1;
						p.y = paths[i][j][m].y - d1 * 0.5;
						p.z = paths[i][j][m].z;
						paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id + 1 && paths[i][j][m].id == paths[i][j][m + 1].id - 1 - x_count) ||
						(paths[i][j][m - 1].id == paths[i][j][m].id - 1 - x_count && paths[i][j][m].id == paths[i][j][m + 1].id + 1)) {//case 20
						p.x = paths[i][j][m].x + d1;
						p.y = paths[i][j][m].y + d1 * 0.5;
						p.z = paths[i][j][m].z;
						paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id + x_count - 1 && paths[i][j][m].id == paths[i][j][m + 1].id - 1 - x_count) ||
						(paths[i][j][m - 1].id == paths[i][j][m].id - 1 - x_count && paths[i][j][m].id == paths[i][j][m + 1].id + x_count - 1)) {//case 21
						p.x = paths[i][j][m].x;
						p.y = paths[i][j][m].y - d1 * sqrt(2) * 0.5;
						p.z = paths[i][j][m].z;
						paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id + x_count + 1 && paths[i][j][m].id == paths[i][j][m + 1].id + 1 - x_count) ||
						(paths[i][j][m - 1].id == paths[i][j][m].id + 1 - x_count && paths[i][j][m].id == paths[i][j][m + 1].id + x_count + 1)) {//case 22
						p.x = paths[i][j][m].x;
						p.y = paths[i][j][m].y + d1 * sqrt(2) * 0.5;
						p.z = paths[i][j][m].z;
						paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id + x_count + 1 && paths[i][j][m].id == paths[i][j][m + 1].id - 1 + x_count) ||
						(paths[i][j][m - 1].id == paths[i][j][m].id - 1 + x_count && paths[i][j][m].id == paths[i][j][m + 1].id + x_count + 1)) {//case 23
						p.x = paths[i][j][m].x - d1 * sqrt(2) * 0.5;
						p.y = paths[i][j][m].y;
						p.z = paths[i][j][m].z;
						paths_buffer1.push_back(p);
						}
						else if ((paths[i][j][m - 1].id == paths[i][j][m].id - x_count + 1 && paths[i][j][m].id == paths[i][j][m + 1].id - 1 - x_count) ||
						(paths[i][j][m - 1].id == paths[i][j][m].id - 1 - x_count && paths[i][j][m].id == paths[i][j][m + 1].id - x_count + 1)) {//case 24
						p.x = paths[i][j][m].x + d1 * sqrt(2) * 0.5;
						p.y = paths[i][j][m].y;
						p.z = paths[i][j][m].z;
						paths_buffer1.push_back(p);
						}
					}
				}
			}
			//paths_buffer1.push_back(paths_buffer1[0]);//终点与起点连接
			paths_buffer2.push_back(paths_buffer1);
			paths_buffer1.clear();
		}
		paths_buffer3.push_back(paths_buffer2);
		paths_buffer2.clear();
	}
	lines = paths;
	paths.clear();
	paths = paths_buffer3;
	paths_buffer3.clear();
	paths_buffer1.clear();
	for (int i = 0; i < paths.size(); i++) {
		for (int j = 0; j < paths[i].size(); j++) {
			for (int m = 0; m < paths[i][j].size() - 1; m += 2) {
				p.x = (paths[i][j][m].x + paths[i][j][m + 1].x) * 0.5;
				p.y = (paths[i][j][m].y + paths[i][j][m + 1].y) * 0.5;
				p.z = (paths[i][j][m].z + paths[i][j][m + 1].z) * 0.5;
				paths_buffer1.push_back(p);
				paths[i][j].insert(paths[i][j].begin() + m + 1, paths_buffer1.begin(), paths_buffer1.end());
				paths_buffer1.clear();

			}
		}
	}
}
void underfill_optimize() {
	path p;
	for (int i = 0; i < model.size(); i++) {
		for (int j = 0; j < model[i].size(); j++) {
			if (model[i][j].b == 0) {
				model[i][j].b++;
				for (int n = 0; n < 4; n++) {
					if (n == 0) {
						p.x = model[i][j].x - d1 * 0.5;
						p.y = model[i][j].y + d1 * 0.5;
						p.z = model[i][j].z;
						paths_buffer1.push_back(p);
					}
					else if (n == 1) {
						p.x = model[i][j].x + d1 * 0.5;
						p.y = model[i][j].y + d1 * 0.5;
						p.z = model[i][j].z;
						paths_buffer1.push_back(p);
					}
					else if (n == 2) {
						p.x = model[i][j].x + d1 * 0.5;
						p.y = model[i][j].y - d1 * 0.5;
						p.z = model[i][j].z;
						paths_buffer1.push_back(p);
					}
					else if (n == 3) {
						p.x = model[i][j].x - d1 * 0.5;
						p.y = model[i][j].y - d1 * 0.5;
						p.z = model[i][j].z;
						paths_buffer1.push_back(p);
					}
				}
				
				paths_buffer1.push_back(paths_buffer1[0]);
				paths[i].push_back(paths_buffer1);
				paths_buffer1.clear();
				p.x = model[i][j].x;
				p.y = model[i][j].y;
				p.z = model[i][j].z;
				paths_buffer1.push_back(p);
				lines[i].push_back(paths_buffer1);
				paths_buffer1.clear();
			}
		}
	}
}

void Continuous_path() {

	for (int i = 0; i < paths.size(); i++) {
		for (int j = 1; j < paths[i].size(); j++) {
			for (int n = 0; n < paths[i][j].size() - 1; n++) {
				for (int m = 0; m < paths[i][0].size() - 1; m++) {
					//上下
					if (paths[i][j][n].x != paths[i][j][n + 1].x && paths[i][j][n].y == paths[i][j][n + 1].y&&paths[i][0][m].x != paths[i][0][m + 1].x && paths[i][0][m].y == paths[i][0][m + 1].y) {
						//顺序
						if (paths[i][j][n].x == paths[i][0][m].x && paths[i][j][n + 1].x == paths[i][0][m + 1].x && abs(paths[i][j][n].y - paths[i][0][m].y) > d1 - 0.01 && abs(paths[i][j][n].y - paths[i][0][m].y) < d1 + 0.01 && abs(paths[i][j][n + 1].y - paths[i][0][m + 1].y) > d1-0.01 && abs(paths[i][j][n + 1].y - paths[i][0][m + 1].y) < d1 + 0.01) {
							reverse(paths[i][j].begin(), paths[i][j].begin() + n + 1);
							reverse(paths[i][j].begin() + n + 1, paths[i][j].end());
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin(), paths[i][j].end());
							paths[i].erase(paths[i].begin() + j);
							j = 0; break;
						}
						////逆序
						else if (paths[i][j][n].x == paths[i][0][m + 1].x && paths[i][j][n+1].x == paths[i][0][m].x && abs(paths[i][j][n].y - paths[i][0][m + 1].y) > d1 - 0.01 && abs(paths[i][j][n].y - paths[i][0][m + 1].y) < d1 + 0.01 && abs(paths[i][j][n+1].y - paths[i][0][m].y) > d1 - 0.01 && abs(paths[i][j][n + 1].y - paths[i][0][m].y) < d1 + 0.01) {
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin(), paths[i][j].begin() + n + 1);
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin() + n + 1, paths[i][j].end());
							paths[i].erase(paths[i].begin() + j);
							j = 0; break;
						}
					}
					//左右
					else if (paths[i][j][n].y != paths[i][j][n + 1].y && paths[i][j][n].x == paths[i][j][n + 1].x&& paths[i][0][m].y != paths[i][0][m + 1].y && paths[i][0][m].x == paths[i][0][m + 1].x) {
						//顺序
						if (paths[i][j][n].y == paths[i][0][m].y && paths[i][j][n + 1].y == paths[i][0][m + 1].y && abs(paths[i][j][n].x - paths[i][0][m].x) > d1 -0.0l && abs(paths[i][j][n].x - paths[i][0][m].x) < d1 + 0.0l && abs(paths[i][j][n + 1].x - paths[i][0][m + 1].x) > d1 - 0.01 && abs(paths[i][j][n + 1].x - paths[i][0][m + 1].x) < d1 + 0.01) {
							reverse(paths[i][j].begin(), paths[i][j].begin() + n + 1);
							reverse(paths[i][j].begin() + n + 1, paths[i][j].end());
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin(), paths[i][j].end());
							paths[i].erase(paths[i].begin() + j);
							j = 0; break;
						}
						//逆序
						else if (paths[i][j][n].y == paths[i][0][m + 1].y && paths[i][j][n+1].y == paths[i][0][m].y && abs(paths[i][j][n].x - paths[i][0][m + 1].x) > d1-0.001&& abs(paths[i][j][n].x - paths[i][0][m + 1].x) < d1 + 0.001 && abs(paths[i][j][n+1].x - paths[i][0][m].x) > d1-0.001&& abs(paths[i][j][n + 1].x - paths[i][0][m].x) < d1 + 0.001) {
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin(), paths[i][j].begin() + n + 1);
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin() + n + 1, paths[i][j].end());
							paths[i].erase(paths[i].begin() + j);
							j = 0; break;
						}
					}
					//斜
					else if (paths[i][j][n].y != paths[i][j][n + 1].y && paths[i][j][n].x != paths[i][j][n + 1].x && paths[i][0][m].y != paths[i][0][m + 1].y && paths[i][0][m].x != paths[i][0][m + 1].x) {
						//顺序
						if (abs(paths[i][j][n].x - paths[i][0][m].x) < (sqrt(2) * d2 * 0.5) + 0.01 && abs(paths[i][j][n].x - paths[i][0][m].x) > (sqrt(2) * d2 * 0.5) - 0.01 && abs(paths[i][j][n].y - paths[i][0][m].y) < (sqrt(2) * d2 * 0.5) + 0.01 && abs(paths[i][j][n].y - paths[i][0][m].y) > (sqrt(2) * d2 * 0.5) - 0.01 &&
							abs(paths[i][j][n+1].x - paths[i][0][m+1].x) < (sqrt(2) * d2 * 0.5) + 0.01 && abs(paths[i][j][n+1].x - paths[i][0][m+1].x) > (sqrt(2) * d2 * 0.5) - 0.01 && abs(paths[i][j][n+1].y - paths[i][0][m+1].y) < (sqrt(2) * d2 * 0.5) + 0.01 && abs(paths[i][j][n+1].y - paths[i][0][m+1].y) > (sqrt(2) * d2 * 0.5) - 0.01) {
							paths[i][0][m + 1].x = (paths[i][0][m].x + paths[i][0][m + 1].x) * 0.5;
							paths[i][0][m + 1].y = (paths[i][0][m].y + paths[i][0][m + 1].y) * 0.5;
							paths[i][j][n + 1].x = (paths[i][j][n].x + paths[i][j][n + 1].x) * 0.5;
							paths[i][j][n + 1].y = (paths[i][j][n].y + paths[i][j][n + 1].y) * 0.5;
							reverse(paths[i][j].begin(), paths[i][j].begin() + n + 1);
							reverse(paths[i][j].begin() + n + 1, paths[i][j].end());
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin(), paths[i][j].end());
							paths[i].erase(paths[i].begin() + j);
							j = 0; break;
						}
						//逆序
						else if (abs(paths[i][j][n].x - paths[i][0][m + 1].x) < (sqrt(2) * d2 * 0.5) + 0.1 && abs(paths[i][j][n].x - paths[i][0][m + 1].x) > (sqrt(2) * d2 * 0.5) - 0.1 && abs(paths[i][j][n].y - paths[i][0][m + 1].y) < (sqrt(2) * d2 * 0.5) + 0.1 && abs(paths[i][j][n].y - paths[i][0][m + 1].y) > (sqrt(2) * d2 * 0.5) - 0.1&&
							abs(paths[i][j][n+1].x - paths[i][0][m].x) < (sqrt(2) * d2 * 0.5) + 0.1 && abs(paths[i][j][n+1].x - paths[i][0][m].x) > (sqrt(2) * d2 * 0.5) - 0.1 && abs(paths[i][j][n+1].y - paths[i][0][m].y) < (sqrt(2) * d2 * 0.5) + 0.1 && abs(paths[i][j][n+1].y - paths[i][0][m].y) > (sqrt(2) * d2 * 0.5) - 0.1) {
							paths[i][0][m + 1].x = (paths[i][0][m].x + paths[i][0][m + 1].x) * 0.5;
							paths[i][0][m + 1].y = (paths[i][0][m].y + paths[i][0][m + 1].y) * 0.5;
							paths[i][j][n].x = (paths[i][j][n].x + paths[i][j][n + 1].x) * 0.5;
							paths[i][j][n].y = (paths[i][j][n].y + paths[i][j][n + 1].y) * 0.5;
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin(), paths[i][j].begin() + n + 1);//左开右闭
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin() + n + 1, paths[i][j].end());
							paths[i].erase(paths[i].begin() + j);
							j = 0; break;
						}
					}
				}
				if (j == 0) break;
			}
		}
	}


	//////孤立点连接
	for (int i = 0; i < paths.size(); i++) {
		for (int j = 1; j < paths[i].size(); j++) {
			for (int n = 0; n < paths[i][j].size() - 1; n++) {
				for (int m = 0; m < paths[i][0].size() - 1; m++) {
					//上下
					if (paths[i][j][n].x != paths[i][j][n + 1].x && paths[i][j][n].y == paths[i][j][n + 1].y && paths[i][0][m].x != paths[i][0][m + 1].x && paths[i][0][m].y == paths[i][0][m + 1].y) {
						//顺序
						if (abs(paths[i][j][n].x - paths[i][0][m].x) < (d1 * 0.5) + 0.4&& abs(paths[i][j][n].x - paths[i][0][m].x) > (d1 * 0.5) - 0.4 && abs(paths[i][j][n + 1].x - paths[i][0][m + 1].x) < (d1*0.5) + 0.4 && abs(paths[i][j][n + 1].x - paths[i][0][m + 1].x) > (d1 * 0.5) - 0.4 && abs(paths[i][j][n].y - paths[i][0][m].y) > d1 - 0.1 && abs(paths[i][j][n].y - paths[i][0][m].y) < d1 + 0.1 && abs(paths[i][j][n + 1].y - paths[i][0][m + 1].y) > d1 - 0.1 && abs(paths[i][j][n + 1].y - paths[i][0][m + 1].y) < d1 + 0.1) {
							reverse(paths[i][j].begin(), paths[i][j].begin() + n + 1);
							reverse(paths[i][j].begin() + n + 1, paths[i][j].end());
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin(), paths[i][j].end());
							paths[i].erase(paths[i].begin() + j);
							j = 0; break;
						}
						////逆序
						else if (abs(paths[i][j][n].x - paths[i][0][m+1].x) < (d1 * 0.5) + 0.4 && abs(paths[i][j][n].x - paths[i][0][m+1].x) > (d1 * 0.5) - 0.4 && abs(paths[i][j][n + 1].x - paths[i][0][m].x) < (d1 * 0.5) + 0.4 && abs(paths[i][j][n + 1].x - paths[i][0][m].x) > (d1 * 0.5) - 0.4 && abs(paths[i][j][n].y - paths[i][0][m + 1].y) > d1 - 0.1 && abs(paths[i][j][n].y - paths[i][0][m + 1].y) < d1 + 0.1 && abs(paths[i][j][n + 1].y - paths[i][0][m].y) > d1 - 0.1 && abs(paths[i][j][n + 1].y - paths[i][0][m].y) < d1 + 0.1) {
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin(), paths[i][j].begin() + n + 1);
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin() + n + 1, paths[i][j].end());
							paths[i].erase(paths[i].begin() + j);
							j = 0; break;
						}
					}
					//左右
					else if (paths[i][j][n].y != paths[i][j][n + 1].y && paths[i][j][n].x == paths[i][j][n + 1].x && paths[i][0][m].y != paths[i][0][m + 1].y && paths[i][0][m].x == paths[i][0][m + 1].x) {
						//顺序
						if (abs(paths[i][j][n].y - paths[i][0][m].y) < (d1 * 0.5) + 0.4 && abs(paths[i][j][n].y - paths[i][0][m].y) > (d1 * 0.5) - 0.4 && abs(paths[i][j][n + 1].y - paths[i][0][m + 1].y) < (d1 * 0.5) + 0.4 && abs(paths[i][j][n + 1].y - paths[i][0][m + 1].y) > (d1 * 0.5) - 0.4 && abs(paths[i][j][n].x - paths[i][0][m].x) > d1 - 0.l && abs(paths[i][j][n].x - paths[i][0][m].x) < d1 + 0.l && abs(paths[i][j][n + 1].x - paths[i][0][m + 1].x) > d1 - 0.1 && abs(paths[i][j][n + 1].x - paths[i][0][m + 1].x) < d1 + 0.1) {
							reverse(paths[i][j].begin(), paths[i][j].begin() + n + 1);
							reverse(paths[i][j].begin() + n + 1, paths[i][j].end());
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin(), paths[i][j].end());
							paths[i].erase(paths[i].begin() + j);
							j = 0; break;
						}
						//逆序
						else if (abs(paths[i][j][n].y - paths[i][0][m + 1].y) < (d1 * 0.5) + 0.4 && abs(paths[i][j][n].y - paths[i][0][m + 1].y) > (d1 * 0.5) - 0.4 && abs(paths[i][j][n + 1].y - paths[i][0][m].y) < (d1 * 0.5) + 0.4 && abs(paths[i][j][n + 1].y - paths[i][0][m].y) > (d1 * 0.5) - 0.4 && abs(paths[i][j][n].x - paths[i][0][m + 1].x) > d1 - 0.1 && abs(paths[i][j][n].x - paths[i][0][m + 1].x) < d1 + 0.1 && abs(paths[i][j][n + 1].x - paths[i][0][m].x) > d1 - 0.1 && abs(paths[i][j][n + 1].x - paths[i][0][m].x) < d1 + 0.1) {
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin(), paths[i][j].begin() + n + 1);
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin() + n + 1, paths[i][j].end());
							paths[i].erase(paths[i].begin() + j);
							j = 0; break;
						}
					}
					//斜
					else if (paths[i][j][n].y != paths[i][j][n + 1].y && paths[i][j][n].x != paths[i][j][n + 1].x && paths[i][0][m].y != paths[i][0][m + 1].y && paths[i][0][m].x != paths[i][0][m + 1].x) {
						//顺序
						if (abs(paths[i][j][n].x - paths[i][0][m].x) < (sqrt(2) * d2 * 0.5) + 0.4 && abs(paths[i][j][n].x - paths[i][0][m].x) > (sqrt(2) * d2 * 0.5) - 0.4 && abs(paths[i][j][n].y - paths[i][0][m].y) < (sqrt(2) * d2 * 0.5) + 0.4 && abs(paths[i][j][n].y - paths[i][0][m].y) > (sqrt(2) * d2 * 0.5) - 0.4 &&
							abs(paths[i][j][n + 1].x - paths[i][0][m + 1].x) < (sqrt(2) * d2 * 0.5) + 0.4 && abs(paths[i][j][n + 1].x - paths[i][0][m + 1].x) > (sqrt(2) * d2 * 0.5) - 0.4 && abs(paths[i][j][n + 1].y - paths[i][0][m + 1].y) < (sqrt(2) * d2 * 0.5) + 0.4 && abs(paths[i][j][n + 1].y - paths[i][0][m + 1].y) > (sqrt(2) * d2 * 0.5) - 0.4) {
							paths[i][0][m + 1].x = (paths[i][0][m].x + paths[i][0][m + 1].x) * 0.5;
							paths[i][0][m + 1].y = (paths[i][0][m].y + paths[i][0][m + 1].y) * 0.5;
							paths[i][j][n + 1].x = (paths[i][j][n].x + paths[i][j][n + 1].x) * 0.5;
							paths[i][j][n + 1].y = (paths[i][j][n].y + paths[i][j][n + 1].y) * 0.5;
							reverse(paths[i][j].begin(), paths[i][j].begin() + n + 1);
							reverse(paths[i][j].begin() + n + 1, paths[i][j].end());
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin(), paths[i][j].end());
							paths[i].erase(paths[i].begin() + j);
							j = 0; break;
						}
						//逆序
						else if (abs(paths[i][j][n].x - paths[i][0][m + 1].x) < (sqrt(2) * d2 * 0.5) + 0.4 && abs(paths[i][j][n].x - paths[i][0][m + 1].x) > (sqrt(2) * d2 * 0.5) - 0.4 && abs(paths[i][j][n].y - paths[i][0][m + 1].y) < (sqrt(2) * d2 * 0.5) + 0.4 && abs(paths[i][j][n].y - paths[i][0][m + 1].y) > (sqrt(2) * d2 * 0.5) - 0.4 &&
							abs(paths[i][j][n + 1].x - paths[i][0][m].x) < (sqrt(2) * d2 * 0.5) + 0.4 && abs(paths[i][j][n + 1].x - paths[i][0][m].x) > (sqrt(2) * d2 * 0.5) - 0.4 && abs(paths[i][j][n + 1].y - paths[i][0][m].y) < (sqrt(2) * d2 * 0.5) + 0.4 && abs(paths[i][j][n + 1].y - paths[i][0][m].y) > (sqrt(2) * d2 * 0.5) - 0.4) {
							paths[i][0][m + 1].x = (paths[i][0][m].x + paths[i][0][m + 1].x) * 0.5;
							paths[i][0][m + 1].y = (paths[i][0][m].y + paths[i][0][m + 1].y) * 0.5;
							paths[i][j][n].x = (paths[i][j][n].x + paths[i][j][n + 1].x) * 0.5;
							paths[i][j][n].y = (paths[i][j][n].y + paths[i][j][n + 1].y) * 0.5;
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin(), paths[i][j].begin() + n + 1);//左开右闭
							paths[i][0].insert(paths[i][0].begin() + m + 1, paths[i][j].begin() + n + 1, paths[i][j].end());
							paths[i].erase(paths[i].begin() + j);
							j = 0; break;
						}
					}
				}
				if (j == 0) break;
			}
		}
	}
	for (int i = 0; i < paths.size(); i++) {
		for (int j = 0; j < paths[i].size(); j++) {
			paths[i][j].push_back(paths[i][j][0]);
		}
	}
	//删除直线的中间坐标点
	for (int i = 0; i < paths.size(); i++) {
		for (int j = 0; j < paths[i].size(); j++) {
			if (paths[i][j].size() > 2) {
				for (int m = 2; m < paths[i][j].size(); m++) {
					if (((paths[i][j][m].x==paths[i][j][m-1].x&&paths[i][j][m-1].x==paths[i][j][m-2].x)||(paths[i][j][m].y == paths[i][j][m - 1].y && paths[i][j][m - 1].y == paths[i][j][m - 2].y))
					) {
						paths[i][j].erase(paths[i][j].begin() + m - 1);
						m = m - 1;
					}
				}
			}

		}
	}
}

//void DFS(int i,int j) {
//	int m;
//	visited[i][j] = 1;//第i层的第j个路径已被访问
//	arr_buffer1.push_back(j);
//	for (m = 0; m < G[i].size(); m++) {
//		if (G[i][j][m] == 1 && visited[i][m] != 1) {
//			DFS(i, m);
//		}
//	}
//}
//
//void DFSTraverse() {//利用深度优先遍历算法打通路径实现一笔画
//	///////////////////生成邻接矩阵
//	int b;
//	for (int i = 0; i < lines.size(); i++) {
//		for (int j = 0; j < lines[i].size(); j++) {
//			for (int m = 0; m < lines[i].size(); m++) {
//				b = 0;//b为是否可连接
//				if (j == m) {
//					G_buffer1.push_back(b);
//					continue;
//				}
//				for (int k = 0; k < lines[i][j].size()-1; k++) {
//					for (int n = 0; n < lines[i][m].size(); n++) {
//						if (abs(lines[i][j][k].id - lines[i][j][k + 1].id) == (1 || x_count) && abs(lines[i][j][k].id - lines[i][j][n].id) == (1 || x_count)) {
//							b = 1;
//							break;
//						}else if (abs(lines[i][j][k].id - lines[i][j][k + 1].id) == x_count - 1 && abs(lines[i][j][k].id - lines[i][j][n].id) == x_count - 1) {
//							b = 1;
//							break;
//						}
//					}
//					if (b == 1)
//						break;
//				}
//				//for (int n = 0; n < lines[i][m].size(); n++) {
//				//	if (abs(lines[i][j][0].id - lines[i][j][1].id) == 1 || abs(lines[i][j][0].id - lines[i][j][1].id) == x_count || abs(lines[i][j][lines[i][j].size() - 1].id - lines[i][j][lines[i][j].size() - 2].id) == 1 || abs(lines[i][j][lines[i][j].size() - 1].id - lines[i][j][lines[i][j].size() - 2].id) == x_count) {
//				//		if (lines[i][j][0].id % x_count == 1 || lines[i][j][lines[i][j].size() - 1].id % x_count == 1) {//左边界
//				//			if (lines[i][j][0].id + 1 == lines[i][m][n].id || lines[i][j][0].id + x_count == lines[i][m][n].id|| lines[i][j][0].id -x_count == lines[i][m][n].id||
//				//				 lines[i][j][lines[i][j].size() - 1].id + 1 == lines[i][m][n].id || lines[i][j][lines[i][j].size() - 1].id + x_count == lines[i][m][n].id || lines[i][j][lines[i][j].size() - 1].id - x_count == lines[i][m][n].id) {
//				//				b = 1;
//				//				break;
//				//			}
//				//		}
//				//		else if (lines[i][j][0].id % x_count == 0 || lines[i][j][lines[i][j].size() - 1].id % x_count == 0) {//右边界
//				//			if ( lines[i][j][0].id - 1 == lines[i][m][n].id ||  lines[i][j][0].id + x_count == lines[i][m][n].id || lines[i][j][0].id - x_count == lines[i][m][n].id ||
//				//				 lines[i][j][lines[i][j].size() - 1].id - 1 == lines[i][m][n].id ||  lines[i][j][lines[i][j].size() - 1].id + 1 == lines[i][m][n].id || lines[i][j][lines[i][j].size() - 1].id - 1 == lines[i][m][n].id) {
//				//				b = 1;
//				//				break;
//				//			}
//				//		}
//				//		else {
//				//			if (lines[i][j][0].id + 1 == lines[i][m][n].id || lines[i][j][0].id + x_count == lines[i][m][n].id || lines[i][j][0].id - 1 == lines[i][m][n].id ||lines[i][j][0].id - x_count == lines[i][m][n].id ||
//				//				lines[i][j][lines[i][j].size() - 1].id + 1 == lines[i][m][n].id || lines[i][j][lines[i][j].size() - 1].id + x_count == lines[i][m][n].id || lines[i][j][lines[i][j].size() - 1].id - 1 == lines[i][m][n].id || lines[i][j][lines[i][j].size() - 1].id - x_count == lines[i][m][n].id) {
//				//				b = 1;
//				//				break;
//				//			}
//				//		}
//				//	}
//				//}
//				G_buffer1.push_back(b);
//			}
//			G_buffer2.push_back(G_buffer1);
//			G_buffer1.clear();
//		}
//		G.push_back(G_buffer2);
//		G_buffer2.clear();
//	}
//	//初始化所有结点为未访问，其中未访问为0，访问过为1；
//	G_buffer1.clear();
//	for (int i = 0; i < G.size(); i++) {
//		for (int j = 0; j < G[i].size(); j++) {
//			G_buffer1.push_back(0);
//		}
//		visited.push_back(G_buffer1);
//		G_buffer1.clear();
//	}
//	///深度优先遍历算法
//	for (int i = 0; i < G.size(); i++) {
//		for (int j = 0; j < G[i].size(); j++) {
//			if (visited[i][j] != 1) {
//				DFS(i, j);
//			}
//		}
//		arr.push_back(arr_buffer1);
//		arr_buffer1.clear();
//	}
//	for (int i = 0; i < arr[0].size(); i++) {
//		printf("%d,", arr[0][i]);
//	}
//}
//
//
//void Continuous_path() {
//	DFSTraverse();//生成可连接路径序列
//	合并路径
//	paths_buffer1.clear();
//	paths_buffer2.clear();
//	paths_buffer3.clear();//初始化缓存路径
//	for (int i = 0; i < arr.size(); i++) {//根据可连接情况重新排列路径顺序
//		for (int j = 0; j < arr[i].size(); j++) {
//			paths_buffer2.push_back(paths[i][arr[i][j]]);
//		}
//		paths_buffer3.push_back(paths_buffer2);
//		paths_buffer2.clear();
//	}
//	paths.clear();
//	paths = paths_buffer3;
//	paths_buffer3.clear();
//	
//	for (int i = 0; i < paths.size(); i++) {
//		for (int j = 0; j < paths[i].size(); j++) {
//			for (int c = 0; c < paths[i].size(); c++) {
//				for (int m = 0; m < paths[i][j].size() - 1; m++) {
//					for (int n = 0; n < paths[i][c].size() - 1; n++) {
//						if (paths[i][j][m].x != paths[i][j][m + 1].x && paths[i][j][m].y == paths[i][j][m + 1].y) {//水平
//							if (paths[i][c][n].x == paths[i][j][m].x && paths[i][c][n].y == paths[i][j][m].y + d1 &&
//								paths[i][c][n + 1].x == paths[i][j][m + 1].x && paths[i][c][n + 1].y == paths[i][j][m + 1].x + d1) {
//								reverse(paths[i][c].begin(), paths[i][c].begin() + n - 1);
//								reverse(paths[i][c].begin() + n, paths[i][c].end());
//								paths[i][j].insert(paths[i][j].begin() + m - 1, paths[i][c].begin(), paths[i][c].end());
//								paths[i].erase(paths[i].begin() + j);
//								j--;
//								c = 0;
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//}

void main(int argc, char** argv) {
	clock_t start, end;
	start = clock();//开始时间
	printf("正在读取文件........\n");
	ReadFile();
	printf("读取成功！\n");
	printf("正常计算最大主应力.......\n");
	FE_Analysis();
	printf("生成路径中........\n");
	StressPath();
	//clipper();
	double_path();
	underfill_optimize();
	Continuous_path();
	printf("路径生成成功！\n");
	GcodePrint();
	//AbaqusData();
	//MatlabData();
	//fangxiangData();
	end = clock();//结束时间

	int s1, s2, s3, s4, s5;
	s1 = 0; s2 = 0; s3 = 0; s4 = 0; s5 = 0;
	for (int i = 0; i < point.size(); i++) {
		for (int j = 0; j < point[0].size(); j++) {
			for (int m = 0; m < point[0][0].size(); m++) {
				if (point[i][j][m].direction > 0 && point[i][j][m].direction < pi / 8) s1++;
				else if (point[i][j][m].direction >= (double)pi / (double)8 && point[i][j][m].direction < (double)pi * (double)((double)3 / (double)8)) s2++;
				else if (point[i][j][m].direction >= (double)pi * (double)((double)3 / (double)8) && point[i][j][m].direction < (double)pi * (double)((double)5 / (double)8)) s3++;
				else if (point[i][j][m].direction >= (double)pi * (double)((double)5 / (double)8) && point[i][j][m].direction < (double)pi * (double)((double)7 / (double)8)) s4++;
				else if (point[i][j][m].direction >= (double)pi * (double)((double)7 / (double)8) && point[i][j][m].direction < (double)pi) s5++;
				//printf("%d %f %f %f %d %f %f %f %f %f\n", point[i][j][m].id, point[i][j][m].x, point[i][j][m].y, point[i][j][m].z, point[i][j][m].b, point[i][j][m].s11, point[i][j][m].s22, point[i][j][m].s12, point[i][j][m].s_max, point[i][j][m].direction);
			}
		}
	}
	printf("%d %d %d %d %d\n", s1, s2, s3, s4,s5);
	int s=0;
	for (int i = 0; i < model.size(); i++) {
		for (int j = 0; j < model[i].size(); j++) {
			if (model[i][j].b == 1) {
				s++;
			}
		}
	}
	cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;  //输出时间（单位：ｓ）
	printf("第1层总共有%d条路径,单元数量数量为%d\n", paths[0].size(),s);

}
