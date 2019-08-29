#pragma once
#include<iostream>
#include<fstream>
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<algorithm>
#include <typeinfo>
#include<regex>  //正则表达式
#include<string.h>  
#include <atlstr.h>  //包含这个头文件就能使用CString
#include<cstring>
#include<map>  //C++的关联容器能实现一键对应多值(multimap)或者一键对应一个值map
#include<time.h>  //统计查询时间
#include<vector>  //可变数组
#include<ctime> 
#include <boost/regex.hpp> //使用boost正则表达式
#include<io.h>
//#include <regex>


//建立结构体

struct fqfiles   //处理fq的结构体
{
	char id[1024] = {0};
	char *seq = NULL;
	char mark[256] = {0};
	char *qscore = NULL;
	struct fqfiles *pNext = NULL;
};
typedef struct fqfiles fq;
//////////////////
struct queue   //处理fa的结构体
{
	char* gene_id = NULL;//存储基因id
	char * gene_seq = NULL;//存储基因序列
	struct queue *pNext = NULL;//指针域，存储下一个节点的地址
};
typedef struct queue Queue;//简化结构体名称

int reversecomplement(std::string source, char **returnresult);  //反向互补函数
int find_terminal_codon(char *strcodon[], const char *gene_seq); //查找终止密码子
std::string getGeneorgnism(std::string name_acquire); //获取fa文件的物种名
Queue *initialization(); //初始化头节点
//实现数据的插入结构体
void enterQueuesort(Queue **phead, const char *id, std::string gene_seq);  //入队即排序
//数据入链表内存
void enterqueue(Queue **phead, const char *gene_id, std::string gene_seq);//数据入链表内存
// 释放内存
Queue * freequeue(Queue *phead);// 释放内存
//显示数据
void printQueue(Queue *phead); //显示数据
Queue* bubblesort(Queue *phead, char ch);  //冒泡排序，根据基因长度排序

//获取节点数（也就是基因个数）
//写一个实验室载体记录报表软件，用于实验室的载体，抗体以及试剂的分类记录
int getnodenum(Queue *phead); //获取节点数（也就是基因个数）
//根据序列特定片段查找，比如具有某个motif的序列全部找出来（实现） //形式参数是容器的嵌套
void searchGeneMotif(Queue *phead, const char *motif, std::map<std::string, std::vector<std::string>> *myresult); //myresult是传出参数
//查找节点（查找基因）
Queue *searchGeneID(Queue *phead, const char *key); //根据基因名找
//修改链表
void modify(Queue *phead, const char *key, const char *value); //修改链表
//字符串添加，字符串删除，字符查找
//判断两个序列是否相等，或者基因id是否重复，使用strcmp（string.h）C语言的头文件
void deleteDuplicates(Queue *phead); //删除重复的记录
int allseqstatistic(Queue *phead, int filterarrange[]);  //过滤掉特定长度的序列

//综合fa文件数据读取处理
int initial_fa_file(const char *path, Queue **phead);
int mystrlen(const char *str); //自己实现的统计字符串长度，返回值是Int

//////////////////////////////////////////////////////////fq文件处理函数
void fq2fa(const char *path); // 读取fq文件,转换为f\

//check fq seq and qscore length,if wrong would be printed and fixed
void checkfqfile(const char *path);
void checkfqdepth(const char *path); //获取fq相同序列的数量
void randomgrapfqseq(const char *path, unsigned long numofread); //随机获取fq片段
unsigned long getfqnum(const char *path); //获取fq的条数
void get_diff_pos_fq_seq(const char *path, std::vector<unsigned long> &mark); //读取mark数量的read到另外一个文件
void RandomFuncGeneratedNonRepeat(unsigned long originNum, unsigned long needrandnum, std::map<unsigned long, unsigned long> &rand); //获取不重复随机数函数
inline std::vector<std::string> split(std::string str,std::string pattern);  //字符串切割函数，返回一个vector数组
void gff3statistic(const char *pathin, const char *pathout); //
void gff3andgenome(const char *pathingff3, const char *pathingenome, const char *pathout);//获取最长转录本
//mRNA正则表达式
inline bool is_regex(const std::string& object, std::vector<std::string> &myresult, const std::string &str);
//获取转录本信息
int checkgff3(const char *pathingff3);
//生成新的gff3文件
int checkgff3mark(const char *pathingff3, const char *outgff3);
//string 字符串转换为数字后排序
std::vector<int>sortstring(std::vector<std::string> &in);
//
int checkgff3update(const char *pathingff3, const char *pathinfa);
//latest version
int checkgff3mark_update(const char *pathingff3, const char *outgff3);

//获取指定路径下的所有文件
int getAllFiles(std::string path, std::vector<std::string>& files);
