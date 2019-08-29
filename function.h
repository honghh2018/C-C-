#pragma once
#include<iostream>
#include<fstream>
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<algorithm>
#include <typeinfo>
#include<regex>  //������ʽ
#include<string.h>  
#include <atlstr.h>  //�������ͷ�ļ�����ʹ��CString
#include<cstring>
#include<map>  //C++�Ĺ���������ʵ��һ����Ӧ��ֵ(multimap)����һ����Ӧһ��ֵmap
#include<time.h>  //ͳ�Ʋ�ѯʱ��
#include<vector>  //�ɱ�����
#include<ctime> 
#include <boost/regex.hpp> //ʹ��boost������ʽ
#include<io.h>
//#include <regex>


//�����ṹ��

struct fqfiles   //����fq�Ľṹ��
{
	char id[1024] = {0};
	char *seq = NULL;
	char mark[256] = {0};
	char *qscore = NULL;
	struct fqfiles *pNext = NULL;
};
typedef struct fqfiles fq;
//////////////////
struct queue   //����fa�Ľṹ��
{
	char* gene_id = NULL;//�洢����id
	char * gene_seq = NULL;//�洢��������
	struct queue *pNext = NULL;//ָ���򣬴洢��һ���ڵ�ĵ�ַ
};
typedef struct queue Queue;//�򻯽ṹ������

int reversecomplement(std::string source, char **returnresult);  //���򻥲�����
int find_terminal_codon(char *strcodon[], const char *gene_seq); //������ֹ������
std::string getGeneorgnism(std::string name_acquire); //��ȡfa�ļ���������
Queue *initialization(); //��ʼ��ͷ�ڵ�
//ʵ�����ݵĲ���ṹ��
void enterQueuesort(Queue **phead, const char *id, std::string gene_seq);  //��Ӽ�����
//�����������ڴ�
void enterqueue(Queue **phead, const char *gene_id, std::string gene_seq);//�����������ڴ�
// �ͷ��ڴ�
Queue * freequeue(Queue *phead);// �ͷ��ڴ�
//��ʾ����
void printQueue(Queue *phead); //��ʾ����
Queue* bubblesort(Queue *phead, char ch);  //ð�����򣬸��ݻ��򳤶�����

//��ȡ�ڵ�����Ҳ���ǻ��������
//дһ��ʵ���������¼�������������ʵ���ҵ����壬�����Լ��Լ��ķ����¼
int getnodenum(Queue *phead); //��ȡ�ڵ�����Ҳ���ǻ��������
//���������ض�Ƭ�β��ң��������ĳ��motif������ȫ���ҳ�����ʵ�֣� //��ʽ������������Ƕ��
void searchGeneMotif(Queue *phead, const char *motif, std::map<std::string, std::vector<std::string>> *myresult); //myresult�Ǵ�������
//���ҽڵ㣨���һ���
Queue *searchGeneID(Queue *phead, const char *key); //���ݻ�������
//�޸�����
void modify(Queue *phead, const char *key, const char *value); //�޸�����
//�ַ�����ӣ��ַ���ɾ�����ַ�����
//�ж����������Ƿ���ȣ����߻���id�Ƿ��ظ���ʹ��strcmp��string.h��C���Ե�ͷ�ļ�
void deleteDuplicates(Queue *phead); //ɾ���ظ��ļ�¼
int allseqstatistic(Queue *phead, int filterarrange[]);  //���˵��ض����ȵ�����

//�ۺ�fa�ļ����ݶ�ȡ����
int initial_fa_file(const char *path, Queue **phead);
int mystrlen(const char *str); //�Լ�ʵ�ֵ�ͳ���ַ������ȣ�����ֵ��Int

//////////////////////////////////////////////////////////fq�ļ�������
void fq2fa(const char *path); // ��ȡfq�ļ�,ת��Ϊf\

//check fq seq and qscore length,if wrong would be printed and fixed
void checkfqfile(const char *path);
void checkfqdepth(const char *path); //��ȡfq��ͬ���е�����
void randomgrapfqseq(const char *path, unsigned long numofread); //�����ȡfqƬ��
unsigned long getfqnum(const char *path); //��ȡfq������
void get_diff_pos_fq_seq(const char *path, std::vector<unsigned long> &mark); //��ȡmark������read������һ���ļ�
void RandomFuncGeneratedNonRepeat(unsigned long originNum, unsigned long needrandnum, std::map<unsigned long, unsigned long> &rand); //��ȡ���ظ����������
inline std::vector<std::string> split(std::string str,std::string pattern);  //�ַ����и��������һ��vector����
void gff3statistic(const char *pathin, const char *pathout); //
void gff3andgenome(const char *pathingff3, const char *pathingenome, const char *pathout);//��ȡ�ת¼��
//mRNA������ʽ
inline bool is_regex(const std::string& object, std::vector<std::string> &myresult, const std::string &str);
//��ȡת¼����Ϣ
int checkgff3(const char *pathingff3);
//�����µ�gff3�ļ�
int checkgff3mark(const char *pathingff3, const char *outgff3);
//string �ַ���ת��Ϊ���ֺ�����
std::vector<int>sortstring(std::vector<std::string> &in);
//
int checkgff3update(const char *pathingff3, const char *pathinfa);
//latest version
int checkgff3mark_update(const char *pathingff3, const char *outgff3);

//��ȡָ��·���µ������ļ�
int getAllFiles(std::string path, std::vector<std::string>& files);
