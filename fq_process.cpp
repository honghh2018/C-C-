#define _CRT_SECURE_NO_WARNINGS
#include "function.h"


void RandomFuncGeneratedNonRepeat(unsigned long originNum, unsigned long needrandnum,std::map<unsigned long,unsigned long> &rand)
{
	using namespace std;
	srand((unsigned int)time(NULL)); //产生随机数，随时间不一样
	vector<unsigned long> temp;
	for (unsigned long i = 0; i < originNum; ++i) {
		temp.push_back(i + 1);
	}
	random_shuffle(temp.begin(), temp.end());  //随机排列
	vector <unsigned long> randnum;
	for (unsigned long i = 0; i < needrandnum; i++)
	{
		randnum.push_back(temp[i]);
	}
	sort(randnum.begin(), randnum.end(), less<unsigned long>()); //排序,升序
	for (unsigned long i = 0; i < randnum.size(); i++)
	{
		rand[randnum[i]] = i;
	}
}


void get_diff_pos_fq_seq(const char *path, std::vector<unsigned long> &mark) //读取给出mark这么多的read,给出起始位置的从起始位置开始读取，没有就从零开始读取
{
	using namespace std;
	printf("get_diff_pos_fq_seq=%d\n", mark.size());
	unsigned long count = 0;
	if (path == nullptr)
	{
		printf("Your input file wasn't exists.\n");
		return;
	}
	const char *outpath = "C:\\Users\\xiaohui\\Desktop\\fq_test1.fa";
	ifstream in(path,ios::in);
	ofstream out(outpath, ios::out);
	string str1, str2, str3, str4;
	if (mark.size()==1)
	{
		while (getline(in, str1))
		{
			//printf("第一个条件\n");
			getline(in, str2);
			getline(in, str3);
			getline(in, str4);
			count++;
			if (count < mark[0])
			{
				out << str1 << "\n" << str2 << "\n" << str3 << "\n" << str4 << endl;
			}
			else if (count ==mark[0])
			{
				out << str1 << "\n" << str2 << "\n" << str3 << "\n" << str4;
				break;
			}
			
		}
	}
	else
	{
		//printf("we in\n");
		while (getline(in, str1))
		{
			getline(in, str2);
			getline(in, str3);
			getline(in, str4);
			count++;
			if (count > mark[0] &&count < mark[1])
			{
				out << str1 << "\n" << str2 << "\n" << str3 << "\n" << str4 << endl;
			}
			else if (count == mark[1])
			{
				out << str1 << "\n" << str2 << "\n" << str3 << "\n" << str4;
				break;
			}
			
		}
	}
	in.close();
	out.close();
}



void fq2fa(const char *path)
{
	using namespace std;
	const char *pathout = "C:\\Users\\xiaohui\\Desktop\\fq_test.fa";
	ifstream in(path, ios::in);
	ofstream out(pathout, ios::out);
	if (!in.is_open() ||!out.is_open())
	{
		printf("open File failured...\n");
		return;
	}

	string line1,line2,line3,line4;
	while (getline(in, line1))
	{
		getline(in, line2);
		getline(in, line3);
		getline(in, line4);  //一次读取四行
		//cout << "id=" << line1 << "\nseq=" << line2 << "\nmark=" << line3 << "\nqscore=" << line4 << endl;
		//cout << "id_length:" << line1.length() << endl;
		//将@标识符替换为>
		line1 = line1.replace(line1.find("@"), 1, ">");
		out << line1 << "\n" << line2 << endl;
	}
	in.close();
	out.close();
}


unsigned long getfqnum(const char *path)  //范围：0~4294967295   0 ~ (2^32 -1)
{
	using namespace std;
	unsigned long count = 0; //记录条数
	if (path == nullptr)
	{
		printf("file inputting wrong...\n");
		return 0;
	}
	ifstream in(path, ios::in);
	if (!in.is_open())
	{
		printf("文件打开失败...\n");
		return 0;
	}
	string line1, line2, line3, line4;
	while (getline(in, line1))
	{
		//cout << line1 << endl;
		getline(in, line2);
		getline(in, line3);
		getline(in, line4);
		count++;
	}
	in.close();
	return count;
}



//check fq seq and qscore length,if wrong would be printed and fixed

void checkfqfile(const char *path)
{
	using namespace std;
	const char *pathlog = "C:\\Users\\xiaohui\\Desktop\\fq_test.log";
	int count = 0; //记录有多少条read被过滤掉
	int count1 = 0; //记录有多少条read
	
	ifstream in(path, ios::in);
	
	if (path == NULL ||!in.is_open())
	{
		return;
	}
	string line1, line2, line3, line4;
	while (getline(in, line1))
	{
		getline(in, line2);
		getline(in, line3);
		getline(in, line4);
		if (line2.length() != line4.length())
		{
			count++;
			continue;
		}
		count1++;
	}
	ofstream log(pathlog, ios::out);
	log << "Total read numbers:" << count1 << endl;
	log << "Filter wrong sequence number:" << count << endl;
	
	if (count >= 1) 
	{
		const char *pathout = "C:\\Users\\xiaohui\\Desktop\\fq_test.fq";

		ofstream out(pathout, ios::out);
		log << "New fq files was generated from fq_test.fq" << endl;
		in.clear(ios::goodbit); //文件指针重置
		in.seekg(ios::beg);
		while (getline(in, line1))
		{
			getline(in, line2);
			getline(in, line3);
			getline(in, line4);
			if (line2.length() != line4.length())
			{
				continue; //跳过这次循环
			}
			out << line1 << "\n" << line2 << "\n" << line3 << "\n" << line4 << endl;
		}
		out.close();
	}
	else
	{
		log << "Your original fq file was correction...\n" << endl;
	}
	log.close();
	in.close();
}

//多少序列是相同的
void checkfqdepth(const char *path)
{
	using namespace std;
	map<string,int> myfq;
	int count = 0;
	ifstream in(path, ios::in);
	if (!in.is_open())
	{
		printf("文件打开失败...\n");
		return;
	}
	string line1, line2, line3, line4;
	while (getline(in, line1))
	{
		getline(in, line2);
		getline(in, line3);
		getline(in, line4);
		myfq[line2]++; //统计相同序列出现的次数
	}
	map<string, int>::iterator it = myfq.begin();
	for(; it != myfq.end(); it++)
	{
		cout << it->first << " " << it->second << endl;
	}

}



//随机取出fq文件的序列条数
void randomgrapfqseq(const char *path, unsigned long numofread)
{
	
	int count = 0; //计数器
	using namespace std;
	if (numofread > getfqnum(path))
	{
		printf("Your numbers inputed out of range,please try a little less...\n");
		return;
	}
	
	std::map<std::string, unsigned long> myfqid;
	//std::vector<std::string> myfqid;
	//std::vector<int> randnum;
	std::map<unsigned long, unsigned long> randnum;
	string line1, line2, line3, line4;
	ifstream in(path, ios::in);
	if (!in.is_open())
	{
		printf("文件打开失败...\n");
		return;
	}
	while (getline(in, line1))
	{
		getline(in, line2);
		getline(in, line3);
		getline(in, line4);
		//myfqid.push_back(line1);
		
		myfqid[line1] = count;
		count++;
	}
	in.clear(ios::goodbit);//在seekg之前需要调用流对象的clear方法，把流的标记清除掉，清除以后就可以正常调用seekg方法了。
	in.seekg(ios::beg); //begin的缩写,将文件指针重新回到头部 //seekp设置输出流文件指针位置，seekg设置输入流文件指针位置
	
	RandomFuncGeneratedNonRepeat(getfqnum(path), numofread, randnum);
	for (std::map<unsigned long, unsigned long>::iterator begin = randnum.begin(); begin != randnum.end(); begin++)
	{
		cout << "key: " << begin->first << "\tvalue: " << begin->second << endl;
	}

	string str1, str2, str3, str4;
	const char *pathout = "C:\\Users\\xiaohui\\Desktop\\randfq.fq";
	ofstream out(pathout, ios::out);
	if (!out.is_open())
	{
		printf("文件打开失败...\n");
		return;
	}
	unsigned long count1 = 1;  //输出文件末尾弃掉换行
	while (getline(in, str1))
	{
		getline(in, str2);
		getline(in, str3);
		getline(in, str4);
		//for (unsigned int i = 0; i < randnum.size(); i++)  //map 的size函数返回map容器的元素个数，也就是键个数
		//{
			if (randnum.find(myfqid[str1]) != randnum.end())  //正确
			{
				count1++;
				//printf("we in\n");
				if (count1 == numofread)
				{
					//printf("we in\n");
					out << str1 << "\n" << str2 << "\n" << str3 << "\n" << str4;
					break;
				}
				else
				{
					out << str1 << "\n" << str2 << "\n" << str3 << "\n" << str4 << endl;
				}
			}
			else
			{
				continue;
			}

		//}
		
	}
	printf("count1=%ld,numofread=%ld\n",count1, numofread);
	in.close();
	out.close();
}









