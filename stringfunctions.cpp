#define _CRT_SECURE_NO_WARNINGS
#include "function.h"
//vc包含的目录：$(VC_IncludePath);$(WindowsSDK_IncludePath);
int mystrlen(const char *str)  //自己实现的字符串长度统计函数
{
	int length = 0;
	if (str == NULL)
	{
		return 0;
	}
	else
	{
		while (*str != '\0')
		{
			length++;
			str++;
		}
	}
	return length;
}

int allseqstatistic(Queue *phead,int filterarrange[])  //过滤掉特定长度的序列
{
	int count_filtered = 0; //记录器,记录被过滤的个数
	int count_reservered = 0; //记录剩下的个数
	const char *pathfilter = "C:\\Users\\xiaohui\\Desktop\\filter_result.txt";
	std::ofstream ouputfilter(pathfilter, std::ios::out);
	if (phead == NULL)
	{
		printf("你输入一个空文件，请检查你文件的路径是否正确...\n");
		return 0;
	}
	if (!ouputfilter.is_open())
	{
		std::cout << "Open file failure!\n";
		return 0;
	}
	Queue *bak = phead;
	while (phead != NULL)
	{
		if (mystrlen(phead->gene_seq) > filterarrange[0] && mystrlen(phead->gene_seq) < filterarrange[1])
		{
			count_reservered++; //统计被保留的序列
		}
		else
		{
			count_filtered++; //统计被过滤掉的序列
		}
		phead = phead->pNext;
	}
	ouputfilter << "#Filter Result:Nnumbers Filtered:" << count_filtered << "; Numbers reservered:" << count_reservered << std::endl;
	ouputfilter << "#Filter arrange: " <<"过滤掉序列长度小于"<< filterarrange[0] << "大于" << filterarrange[1] << "序列"<<std::endl;
	while (bak != NULL)
	{
		if (mystrlen(bak->gene_seq) > filterarrange[0] && mystrlen(bak->gene_seq) < filterarrange[1])
		{
			ouputfilter << ">" << bak->gene_id << "\n" << bak->gene_seq << std::endl;
		}
		bak = bak->pNext;
	}
	return 1; //返回1代表执行成功
}


int reversecomplement(std::string source, char **returnresult)
{
	std::reverse(source.begin(), source.end()); //逆转
	//char* ptemp = const_cast<char *> (source.data()); //去掉const属性
	char *ptemp = (char *)malloc(sizeof(char)*(source.length() + 1)); //分配内存空间，calloc是并初始化为0，molloc不会
	//char *ptemp = new char[tempstr.length() + 1];
	char *bak = ptemp;
	if (ptemp == NULL)
	{
		printf("memory distribution failed...\n");
		return std::string::npos;
	}
	strcpy(ptemp, source.c_str());  //将string类型转变为c语言类型
	//printf("ptemp=%s\n", ptemp);
	//printf("ptemp type=%s",typeid(bak).name());
	for (; *ptemp != '\0'; ptemp++)
	{
		if (*ptemp == 'A' || *ptemp == 'a')
		{
			*ptemp = 'T';
			continue;
		}
		if (*ptemp == 'G' || *ptemp == 'g')
		{
			*ptemp = 'C';
			continue;
		}
		if (*ptemp == 'T' || *ptemp == 't')
		{
			*ptemp = 'A';
			continue;
		}
		if (*ptemp == 'C' || *ptemp == 'c')
		{
			*ptemp = 'G';
			continue;
		}
		if (*ptemp == 'N' || *ptemp == 'n')
		{
			*ptemp = 'N';
			continue; //一定要加不然就有访问冲突
		}
	}
	*returnresult = bak; //传出参数
	//free(ptemp); //不能再函数体内部释放，因为如果释放就会访问冲突,在main中再free
	return 0; //返回0代表成功
}



int find_terminal_codon(char *strcodon[], const char *gene_seq)  //查找终止密码子
{
	if (strcodon == NULL)
	{
		return -1; //failure
	}
	for (char * *start = strcodon; start < strcodon + 3; start++)
	{
		std::string str1(gene_seq);
		size_t strlength = str1.length();
		int find_tail_string_pos = str1.rfind(*start); //需要星号
		if (find_tail_string_pos == strlength - 3)
		{
			return 0;
		}
		else
		{
			continue; //没找到就继续
		}
	}
	return -1; //如果都没找到就返回-1
}


std::string getGeneorgnism(std::string name_acquire) //获取fa文件的物种名
{
	std::string obtained_orgnism; //存储获得的字符串
	if (name_acquire.empty()) //如果是0代表有数据,有数据就是假不执行这个返回操作，1代表无数据
	{
		return "-1"; //failure ，gene_id为空
	}


	std::string match_str(":[[:blank:]]([[:alpha:]]+[[:blank:]][[:alpha:]]+)[[:blank:]]"); //正则表达式
	std::regex rule(match_str);
	std::smatch sm; //获取匹配的字符串
	for (auto it = name_acquire.cbegin(); std::regex_search(it, name_acquire.cend(), sm, std::regex(match_str)); it = sm.suffix().first)
	{
		if (sm.size() != 0)
		{
			obtained_orgnism = sm.str();  //获取带有冒号的匹配字符串
			break;
		}

	}

	return obtained_orgnism;
}


//初始化头节点
Queue *initialization()
{
	return NULL;
}

//实现数据的插入结构体
void enterQueuesort(Queue **phead, const char *id, std::string gene_seq)  //入队即排序
{
	Queue *bak_for_return = NULL;
	//无论是不是空的指针我们这里新new一个内存，用于存放要插入
	Queue *pnewelement = (Queue *)malloc(sizeof(Queue)); //建立对象
	//printf("queue=%p",pnewelement);
	pnewelement->gene_id = (char *)calloc(strlen(id) + 1, sizeof(char)); //分别给指针开辟内存,并且初始化为0
	//printf("gene_id=%p\n", pnewelement->gene_id);
	pnewelement->gene_seq = (char *)calloc(gene_seq.length() + 1, sizeof(char));
	//拷贝字符串
	strcpy(pnewelement->gene_id, id);
	strcpy(pnewelement->gene_seq, gene_seq.data());  //用return试试
	pnewelement->pNext = NULL; //下一个指针等于NULL
	if (*phead == NULL) //如果为空，直接加到后面就行(2019-4-23一直报错的原因就是这里的phead没有给星号，phead==NULL是错误的，*phead==NULL才是对的)
	{
		*phead = pnewelement;
		pnewelement->pNext = NULL;
	}
	else  //传入的头节点已经有数据了
	{
		//从大到小的排序
		if (strlen(pnewelement->gene_seq) > strlen((*phead)->gene_seq))
		{
			pnewelement->pNext = *phead;
			*phead = pnewelement;
		}
		else
		{
			//尾部插入
			Queue *ploop = *phead;
			while (ploop->pNext != NULL) //从头部循环到尾部
			{
				ploop = ploop->pNext;
			}
			if (strlen(pnewelement->gene_seq) <= strlen(ploop->gene_seq))
			{
				ploop->pNext = pnewelement;
				pnewelement->pNext = NULL; //实现尾部大小对换
			}
			//实现中部插入
			Queue *p1, *p2;
			p1 = p2 = NULL;
			p1 = *phead;
			while (p1->pNext != NULL)
			{
				p2 = p1->pNext;
				if (strlen(p1->gene_seq) >= strlen(pnewelement->gene_seq) && strlen(p2->gene_seq) < strlen(pnewelement->gene_seq))
				{
					p1->pNext = pnewelement;
					pnewelement->pNext = p2;
					break;
				}
				p1 = p1->pNext;//向前循环
			}
		}
	}
}


//数据入链表内存
void enterqueue(Queue **phead, const char *gene_id, std::string gene_seq)
{

	if (*phead == NULL)
	{
		//Queue *pnewnode = (Queue *)calloc(1, sizeof(Queue));
		//Queue*pnewnode = (Queue *)malloc(sizeof(Queue));
		Queue *pnewnode = new Queue;
		if (pnewnode == NULL)
		{
			printf("内存分配失败\n");
			return;
		}
		//pnewnode->gene_id = (char *)calloc(strlen(gene_id) + 1, sizeof(char));
		//pnewnode->gene_id = (char *)malloc(sizeof(char)*(strlen(gene_id) + 1));
		pnewnode->gene_id = new char[strlen(gene_id) + 1];
		//pnewnode->gene_seq = (char *)calloc(gene_seq.length() + 1, sizeof(char));
		pnewnode->gene_seq = (char *)malloc(sizeof(char)*(gene_seq.length() + 1));
		pnewnode->gene_seq = new char[gene_seq.length() + 1];
		//memset(pnewnode->gene_id, 0x00, sizeof(char)*(strlen(gene_id) + 1));
		//memset(pnewnode->gene_seq, 0x00, sizeof(char)*(gene_seq.length() + 1));
		strcpy(pnewnode->gene_id, gene_id);
		strcpy(pnewnode->gene_seq, gene_seq.data());
		pnewnode->pNext = NULL;
		*phead = pnewnode;
	}
	else
	{
		//尾部链接上
		Queue *ploop = *phead; //save head node
		while (ploop->pNext != NULL) //
		{
			ploop = ploop->pNext; //一直循环
		}
		//Queue *pnewnode = (Queue *)malloc(sizeof(Queue));
		Queue *pnewnode1 = new Queue;
		if (pnewnode1 == NULL)
		{
			printf("内存分配失败\n");
			return;
		}
		//pnewnode->gene_id = (char *)calloc(strlen(gene_id) + 1, sizeof(char));
		//pnewnode->gene_id = (char *)malloc(sizeof(char)*(strlen(gene_id) + 1));
		pnewnode1->gene_id = new char[strlen(gene_id) + 1];
		pnewnode1->gene_seq = new char[gene_seq.length() + 1];
		//pnewnode->gene_seq = (char *)calloc(gene_seq.length() + 1, sizeof(char));
		//pnewnode->gene_seq = (char *)malloc(sizeof(char)*(gene_seq.length() + 1));
		//memset(pnewnode->gene_id, 0x00, sizeof(char)*(strlen(gene_id) + 1));
		//memset(pnewnode->gene_seq, 0x00, sizeof(char)*(gene_seq.length() + 1));
		strcpy(pnewnode1->gene_id, gene_id);
		strcpy(pnewnode1->gene_seq, gene_seq.data());
		pnewnode1->pNext = NULL;
		ploop->pNext = pnewnode1; //链接上
	}
}

//释放内存
Queue * freequeue(Queue *phead)
{
	Queue *p1, *p2; //防止野指针
	p1 = p2 = NULL;
	p1 = phead; //保留头节点

	//printf("p1->pNext=%s\n",p1->pNext->gene_id);
	//printf("p1->pNext->pNext=%s", p1->pNext->gene_seq);
	while (p1->pNext != NULL) //这里为什么p1==NULL了，怎么p1变成了空指针了？ //p1节点不动，p1结构体的pnext连接上其他的节点，然后删除中间的
	{
		p2 = p1->pNext; //保留下一个节点地址
		p1->pNext = p2->pNext; //释放掉p2节点   //这里是p1->pNext=p2->pNext不是p1=p2->pNext,搞了好久，要注意了
		free(p2->gene_id);
		free(p2->gene_seq);
		delete p2; //释放掉p2内存 
		//printf("\n\n\n");
		//printQueue(phead); //查看删除的中间状态

	}
	free(phead->gene_id);
	free(phead->gene_seq);
	free(phead); //释放头节点
	return NULL; //一定要有这个NULL
}


//显示数据
void printQueue(Queue *phead)
{
	using namespace std;
	const char * path = "C:\\Users\\xiaohui\\Desktop\\debug333.list";
	FILE *pout = fopen(path, "w");
	if (phead == NULL || pout == NULL)
	{
		printf("你打印一个空链表...\n");
		return;
	}
	while (phead != NULL)
	{
		fprintf(pout, "GeneID:%s\nGene_seq:%s\n", phead->gene_id, phead->gene_seq);
		cout << "geneID:" << phead->gene_id << "\n" << "GeneSeq:" << phead->gene_seq << "  length=" << strlen(phead->gene_seq) << endl;
		//pout << "geneID:" << phead->gene_id << "\n" << "GeneSeq:" << phead->gene_seq << "  length=" << strlen(phead->gene_seq) << "\n";
		//printf("%p,%p\n",phead,phead->pNext); //打印节点地址
		phead = phead->pNext;
	}
	fclose(pout);
	//return 1; //执行成功返回1
}

Queue* bubblesort(Queue *phead, char ch)  //冒泡排序，根据基因长度排序
{
	if (phead == NULL || phead->pNext == NULL)
	{
		printf("your are in\n");
		printf("linktable empty or just a element!\n");
		return NULL;
	}
	else
	{
		if (ch == '>')   //从大到小排序
		{
			for (Queue *p1 = phead; p1 != NULL; p1 = p1->pNext) //保存头节点，使用备份循环
			{
				//int p1_len = strlen(p1->gene_seq); //不能这样得到长度后比较，是错误的，必须在if语句中进行长度判断，这里留下来做警示
				//printf("p1_len=%d\n", p1_len);
				//printf("p1_seq=%s\n",p1->gene_seq);
				for (Queue *p2 = phead; p2 != NULL; p2 = p2->pNext)
				{

					//int p2_len = strlen(p2->gene_seq);
					//printf("p2_len=%d\n", p2_len);
					if (strlen(p1->gene_seq) > strlen(p2->gene_seq)) //从大到小排序
					{ //只需要拷贝类的数据存储就行，方法是共享的
						Queue temp;
						temp.gene_id = p1->gene_id; //他娘的指针不用拷贝，字符串才用拷贝
						temp.gene_seq = p1->gene_seq;

						p1->gene_id = p2->gene_id;
						p1->gene_seq = p2->gene_seq;
						p2->gene_id = temp.gene_id;
						p2->gene_seq = temp.gene_seq;
						//strcpy(p1->value, p2->value); //错误的，留作警示
						//strcpy(p1->key, p2->key);
						//strcpy(p2->value, temp.value);
						//strcpy(p2->key, temp.key);
					}
				}
			}
			return phead;
		}
		else if (ch == '<')
		{
			for (Queue *p1 = phead; p1 != NULL; p1 = p1->pNext) //保存头节点，使用备份循环
			{
				//int p1_len = strlen(p1->gene_seq); //不能这样得到长度后比较，是错误的，必须在if语句中进行长度判断，这里留下来做警示
				//printf("p1_len=%d\n", p1_len);
				//printf("p1_seq=%s\n",p1->gene_seq);
				for (Queue *p2 = phead; p2 != NULL; p2 = p2->pNext)
				{

					//int p2_len = strlen(p2->gene_seq);
					//printf("p2_len=%d\n", p2_len);
					if (strlen(p1->gene_seq) < strlen(p2->gene_seq)) //从小到大排序
					{ //只需要拷贝类的数据存储就行，方法是共享的
						Queue temp;
						temp.gene_id = p1->gene_id; //他娘的指针不用拷贝，字符串才用拷贝
						temp.gene_seq = p1->gene_seq;

						p1->gene_id = p2->gene_id;
						p1->gene_seq = p2->gene_seq;
						p2->gene_id = temp.gene_id;
						p2->gene_seq = temp.gene_seq;
						//strcpy(p1->value, p2->value); //错误的，留作警示
						//strcpy(p1->key, p2->key);
						//strcpy(p2->value, temp.value);
						//strcpy(p2->key, temp.key);
					}
				}
			}
			return phead;
		}

	}
	//使用冒泡排序法
	return NULL;
}


//获取节点数（也就是基因个数）
//写一个实验室载体记录报表软件，用于实验室的载体，抗体以及试剂的分类记录
int getnodenum(Queue *phead)
{
	int count = 0; //节点计数器
	if (phead == NULL)
	{
		printf("空链表...\n");
		return 0;
	}
	Queue *pbak = phead;
	while (pbak != NULL)
	{
		count++;
		pbak = pbak->pNext;
	}
	return count;  //统计节点个数
}


//根据序列特定片段查找，比如具有某个motif的序列全部找出来（实现） //形式参数是容器的嵌套
void searchGeneMotif(Queue *phead, const char *motif, std::map<std::string, std::vector<std::string>> *myresult) //myresult是传出参数
{  //容器的嵌套
	using namespace std;
	//multimap<string, string> *myresult = new multimap<string, string>;
	//multimap<string, string> myresult;
	string querystring(motif); //将传入进来的motif转变为string类型
	transform(querystring.begin(), querystring.end(), querystring.begin(), ::toupper); //转换为大写
	if (phead == NULL)
	{
		printf("空链表...\n");
	}
	while (phead != NULL)
	{
		vector<std::string> vectemp;
		vector<std::string> tempvector;
		string subject(phead->gene_seq); //将char *指针指向的字符串转换为string
		transform(subject.begin(), subject.end(), subject.begin(), ::toupper); //转换为大写
		int count = 0;  //计数器
		int begin = -1; //记录位置信息
		while ((begin = subject.find(querystring, begin + 1)) != string::npos)  //find函数是返回找到的第一个的位置，从begin位置开始查找，这里是0。
		{

			//insert(map<string,string>::value_type(cat_gene_id.data(),current->gene_seq));
			count++;
			char str[256] = { 0 };
			sprintf(str, "%d", begin + 1); //从1开始打印给用户看
			string temp(str); //将数字转换为字符串
			tempvector.push_back(temp);
			//vectemp.push_back(temp); //将位置信息存储到vect容器中
			begin = begin + querystring.length();  //移过查找到第一个字符串的后面，重新查找
		}
		if (count != 0) //说明该Gene seq找到motif
		{
			char str1[256] = { 0 };
			sprintf(str1, "%d", count);
			string temp1(str1);
			string geneSeq(phead->gene_seq);
			temp1 = "Number:" + temp1;
			//将位置信息用-->箭头连接到一起组成一个字符串
			string catstring;
			for (size_t i = 1; i < tempvector.size(); ++i)
			{
				tempvector[0] += "-->";
				tempvector[0] += tempvector[i];
			}
			tempvector.resize(1); //重新分配内存大小
			vectemp.push_back(tempvector[0]);
			vectemp.push_back(temp1); //将位置信息存储到vect容器中
			vectemp.push_back(geneSeq);
			//outputsame_seqAndID.insert(map<string,string>::value_type(cat_gene_id.data(),current->gene_seq))
			//(*myresult).insert(make_pair(phead->gene_id, vectemp)); //存储到容器中
			(*myresult).insert(map<string, vector<string>>::value_type(phead->gene_id, vectemp));
			//(*myresult).insert(pair<string, string>(phead->gene_id, phead->gene_seq));
			//(*myresult).insert(pair<string, string>(phead->gene_id, temp1)); //在容器中的位置看插入时的顺序
		}
		phead = phead->pNext; //循环
	}
}


//查找节点（查找基因）
Queue *searchGeneID(Queue *phead, const char *key) //根据基因名找
{
	if (phead == NULL)
	{
		printf("空链表...\n");
		return NULL;
	}
	char keyarray[1024] = { 0 };
	strcpy(keyarray, key);
	//keyarray[strlen(key) - 1] = '\0'; //去除\n,这里是错误的，会把传入的字符串弃掉倒数第一个字符
	for (size_t i = 0; i < strlen(keyarray); i++)  //小写转大写
	{
		keyarray[i] = toupper(keyarray[i]); //26个英文字母会转换，其他字符不会
	}
	printf("searchkey=%s\n", keyarray);
	while (phead != NULL)
	{
		char pheadarray[1024] = { 0 };
		strcpy(pheadarray, phead->gene_id);
		printf("pheadarray=%s\n", pheadarray);
		for (size_t i = 0; i < strlen(pheadarray); i++)  //小写转大写,因为strlen(p1)返回Unsigned int,而你定义的i也必须是unsigned否者报警告，这里使用size_t
		{
			pheadarray[i] = toupper(pheadarray[i]); //26个英文字母会转换，其他字符不会
		}
		if (strcmp(keyarray, pheadarray) == 0) //两个字符串相等
		{
			return phead; //返回当前指针节点的地址
		}
		phead = phead->pNext; //循环遍历
		//delete[]p1;  //找到返回了就不能释放堆内存了
	}
	return NULL; //如果查找完了还没有找到返回NULL
}


//修改链表
void modify(Queue *phead, const char *key, const char *value)
{
	if (phead == NULL)
	{
		return;
	}
	Queue *psearch = searchGeneID(phead, key);
	if (psearch == NULL)
	{
		printf("没找到..\n");
	}
	else
	{
		strcpy(psearch->gene_id, key);
		printf("修改成功..\n"); //这个函数还没有运行过
	}
}


//字符串添加，字符串删除，字符查找
//判断两个序列是否相等，或者基因id是否重复，使用strcmp（string.h）C语言的头文件

void deleteDuplicates(Queue *phead)
{
	using namespace std;
	//map元素默认按照键的升序排序
	map<string, string> outputsame_seqAndID; //用来存储输出字符串的键和值（gene_id和gene_seq）他们是键不同值相同
	int count = 0; //记录删除的重复数量
	int count1 = 0; //记录id不同序列相同的个数
	if (phead == NULL)
	{
		return;
	}
	else
	{
		Queue *current = phead;
		Queue *nextnode = NULL;
		while (current != NULL)
		{
			nextnode = current->pNext;
			if (nextnode == NULL)
			{
				break; //循环到结尾就跳出循环
			}
			if (strcmp(current->gene_id, nextnode->gene_id) != 0 && nextnode && strcmp(current->gene_seq, nextnode->gene_seq) == 0) //输出序列一样但是gene ID不一样的
			{
				string cat_gene_id(current->gene_id);
				cat_gene_id += "------>";
				cat_gene_id += nextnode->gene_id;
				outputsame_seqAndID.insert(map<string, string>::value_type(cat_gene_id.data(), current->gene_seq));
				count1++;
			}
			if (strlen(current->gene_seq) == strlen(nextnode->gene_seq) && strcmp(current->gene_id, nextnode->gene_id) == 0 && nextnode && strcmp(current->gene_seq, nextnode->gene_seq) == 0) //strcpm如果相等字符串返回0
			{
				count++;
				current->pNext = nextnode->pNext;
				//printf("current_seq=%s,nextnode=%s\n", current->gene_seq, nextnode->gene_seq);
				free(nextnode);
			}
			else
			{
				current = nextnode;  //循环
			}
		}

	}
	if (!outputsame_seqAndID.empty()) //成立了才创建文件，防止出现垃圾文件,empty函数，非空返回0，空返回1，这里取反
	{
		const char *path1 = "C:\\Users\\xiaohui\\Desktop\\test_equal.fa"; //打开文件句柄输出到该目录
		ofstream outpf(path1, ios::out);
		outpf << "#Same sequences statistic:\n";
		outpf << "#Same sequence numbers:" << count1 << "\tTotal Identical:" << count << "\tDelete Duplicate number:" << count << endl;
		for (map<string, string>::iterator _put2file = outputsame_seqAndID.begin(); _put2file != outputsame_seqAndID.end(); _put2file++)
		{
			outpf << ">" << _put2file->first << "\n" << _put2file->second << endl;
		}

		outpf.close();
	}
}


int initial_fa_file(const char *path,Queue ** phead)
{
	using namespace std;
	//综合数据前处理
	fstream pf; //有三个类，ifstream,ofstream和fstream
	pf.open(path, ios::in); //按照输入的方式读取
	//输出流
	fstream pfo;
	const char *outfile = "C:\\Users\\xiaohui\\Desktop\\outputfile2.txt";
	const char *outrevcom = "C:\\Users\\xiaohui\\Desktop\\revcomp.fa"; //输出反向互补序列
	ofstream out(outrevcom, ios::out);
	//ofstream out(outfile, ios::out | ios::app); //普通写入打开或者追加打开//正确的写法
	pfo.open(outfile, ios::out); //将内存数据写入磁盘文件 //也正确
	if (!pf.is_open() || !pfo.is_open() || !out.is_open())  //如果是假打开文件失败
	{
		cout << "Open file failure!\n";
		return -1;
	}

	//获取fa文件的物种名（目前限定是NCBI的fa文件）
	string getname;
	getline(pf, getname, '>'); //读取第一行的gene_id,为空,丢弃
	getline(pf, getname, '>'); //读取第二行
	size_t position_n = getname.find('\n'); //查找第一次换行符出现的位置，为了得到基因的ID;size_t是unsigned int类型
	getname.assign(getname, 0, position_n); //才是对的
	string orgnism_name_obtained = getGeneorgnism(getname);
	pf.clear(ios::goodbit); //文件指针重置
	pf.seekg(ios::beg);
	//输出格式设置,标题
	pfo << "#Statistic File\t" << "Orgnism" << orgnism_name_obtained << endl; //匹配的字符串中得到冒号和空格，这里就不需要冒号和空格
	int gene_num = 0;
	//统计基因个数
	while (!pf.eof()) //没到文件末尾就继续
	{
		char c;
		pf >> c;
		if (c == '>')
		{
			gene_num++;
		}
	}
	pf.clear(ios::goodbit);//在seekg之前需要调用流对象的clear方法，把流的标记清除掉，清除以后就可以正常调用seekg方法了。
	pf.seekg(ios::beg); //begin的缩写,将文件指针重新回到头部 //seekp设置输出流文件指针位置，seekg设置输入流文件指针位置
	/*ios::beg：文件流的起始位置
	ios::cur：文件流的当前位置
	ios::end：文件流的结束位置*/
	//输出文件格式
	pfo << "#Gene_number: " << gene_num << endl;
	pfo << "#Gene_ID\t" << "Gene_length(bp)\t" << "A(bp)" << "\t" << "T(bp)" << "\t" << "G(bp)" << "\t" << "C(bp)" << "\t" << "N(bp)\t" << "GC(%)" << endl;


	string line;
	while (getline(pf, line, '>')) //按照>分隔符号读取文件
	{
		string gene_id; //获取基因id
		string gene_seq; //获取基因序列
		int A, T, G, C, N; //统计碱基数
		A = T = G = C = N = 0;
		double GC_content; //统计GC含量
		GC_content = 0.0;
		//if(line.find('\n')==string::npos);//npos就是-1
		if (line.length() == 0 || line.find(0x0A) == 0)   //用来查找上一行是否是空或者回车,换行符都在0位置，0x0A=‘\n’
		{
			continue; //跳过开头的空行
		}
		size_t pos_n = line.find('\n'); //查找第一次换行符出现的位置，为了得到基因的ID;size_t是unsigned int类型
		//gene_id = line.assign(line, 0, pos_n); 错误的
		gene_id.assign(line, 0, pos_n); //才是对的
		if (size_t pos_space = gene_id.find(0x20))  //找到space 空格
		{
			gene_id.assign(gene_id, 0, pos_space); //截取
		}
		if (size_t pos_tab = gene_id.find(0x09)) //找到水平制表符
		{
			gene_id.assign(gene_id, 0, pos_tab);
		}
		//gene_seq = line.substr(pos_n+1); //获取除了基因名称的剩下序列  //这个可以用
		gene_seq.assign(line, pos_n + 1, line.length() - pos_n); //+1是为了去掉结尾的换行符号  //这个也可以用

		//remove函数在头文件#include<algorithm>中，不包含这个头文件remove函数就不能使用
		//gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), '\n'), gene_seq.end()); 
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x0A), gene_seq.end()); //删除换行符，同上
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x09), gene_seq.end()); //删除“水平制表符”
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x08), gene_seq.end()); //删除backspace
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x20), gene_seq.end()); //删除空格
		int length_seq = gene_seq.length(); //计算序列的长度,并计算是否为3的倍数

		//查看基因序列是否ATG开头  //transform是algorithm头文件中的库
		transform(gene_seq.begin(), gene_seq.end(), gene_seq.begin(), ::toupper); //将字符串全部转换为大写 ，::tolower是小写
		int atg_head = gene_seq.find("ATG"); //起始密码子 //size_t是unsigned int
		//终止密码子UAA、UGA、UAG
		char *termination_codon[3] = { "TAA","TGA","TAG" }; //指针数组，每一个元素指针指向一个字符串(头指针是一个常量地址代表第一个元素)

		int is_find_tail = find_terminal_codon(termination_codon, gene_seq.data());
		//反向互补序列
		char *revp = NULL;
		if (reversecomplement(gene_seq, &revp) != string::npos) //返回-1代表失败
		{
			string revcom(revp); //将传回来的C风格数组转换为C++风格数组
			string rev_gene = gene_id + "_ReverseComplementSeq";
			out << ">" << rev_gene << "\n" << revcom << "\n";  //打印到输出文件
			free(revp); //释放调用函数revcom开的堆内存
		}

		//内存排序之入队即排序
		//enterQueuesort(&phead, gene_id.data(), gene_seq); //接收返回来的堆内存首地址(有问题)
		enterqueue(phead, gene_id.data(), gene_seq);
		//显示数据
		//printQueue(phead);
		//统计字符
		for (string::iterator begin = gene_seq.begin(); begin != gene_seq.end(); begin++) //使用迭代器统计字符数
		{
			if (*begin == 'A' || *begin == 'a')
			{
				A++;
			}
			else if (*begin == 'G' || *begin == 'g')
			{
				G++;
			}
			else if (*begin == 'C' || *begin == 'c')
			{
				C++;
			}
			else if (*begin == 'T' || *begin == 't')
			{
				T++;
			}
			else
			{
				N++;
			}
		}
		string::iterator begin_find = gene_seq.begin();
		string::iterator end_find = gene_seq.end();  //这两个迭代器用于查找字符串

		//char *pbuf = (char *)calloc(gene_seq.length(),sizeof(char));//很重要，这里不要删除
		//strcpy(buf, gene_seq.data());  //data函数获取string的c风格字符串数组 还可以这样转换，strcpy(buf,gene_seq.c_str);c_str是c字符串string的简写
		//strcpy(pbuf, gene_seq.c_str());
		//free(pbuf); //分配和释放内存太快出现问题
		//fprintf(stderr,"c_str=%s\n",buf); //打印到错误输出上
		char gc_str[256] = { 0 };
		GC_content = ((G + C)*1.0 / (A + T + G + C)*1.0)*100.0;
		sprintf(gc_str, "%0.2f", GC_content);

		if (atg_head == 0 && length_seq % 3 == 0) //如果头部有ATG，并且是能被3除尽的就加一个星号
		{
			gene_id += '*'; //在尾部添加一个字符*
		}
		if (is_find_tail == 0) //如果尾部还有终止密码子的话，就用两个星标注
		{
			gene_id += '*';
		}
		pfo << gene_id << "\t" << length_seq << "\t" << A << "\t" << T << "\t" << G << "\t" << C << "\t" << N << "\t" << gc_str << endl;

	}
	pf.close();
	pfo.close();
	out.close(); //关闭文件流
	return 1;
}