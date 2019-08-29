#define _CRT_SECURE_NO_WARNINGS
#include "function.h"
//vc������Ŀ¼��$(VC_IncludePath);$(WindowsSDK_IncludePath);
int mystrlen(const char *str)  //�Լ�ʵ�ֵ��ַ�������ͳ�ƺ���
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

int allseqstatistic(Queue *phead,int filterarrange[])  //���˵��ض����ȵ�����
{
	int count_filtered = 0; //��¼��,��¼�����˵ĸ���
	int count_reservered = 0; //��¼ʣ�µĸ���
	const char *pathfilter = "C:\\Users\\xiaohui\\Desktop\\filter_result.txt";
	std::ofstream ouputfilter(pathfilter, std::ios::out);
	if (phead == NULL)
	{
		printf("������һ�����ļ����������ļ���·���Ƿ���ȷ...\n");
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
			count_reservered++; //ͳ�Ʊ�����������
		}
		else
		{
			count_filtered++; //ͳ�Ʊ����˵�������
		}
		phead = phead->pNext;
	}
	ouputfilter << "#Filter Result:Nnumbers Filtered:" << count_filtered << "; Numbers reservered:" << count_reservered << std::endl;
	ouputfilter << "#Filter arrange: " <<"���˵����г���С��"<< filterarrange[0] << "����" << filterarrange[1] << "����"<<std::endl;
	while (bak != NULL)
	{
		if (mystrlen(bak->gene_seq) > filterarrange[0] && mystrlen(bak->gene_seq) < filterarrange[1])
		{
			ouputfilter << ">" << bak->gene_id << "\n" << bak->gene_seq << std::endl;
		}
		bak = bak->pNext;
	}
	return 1; //����1����ִ�гɹ�
}


int reversecomplement(std::string source, char **returnresult)
{
	std::reverse(source.begin(), source.end()); //��ת
	//char* ptemp = const_cast<char *> (source.data()); //ȥ��const����
	char *ptemp = (char *)malloc(sizeof(char)*(source.length() + 1)); //�����ڴ�ռ䣬calloc�ǲ���ʼ��Ϊ0��molloc����
	//char *ptemp = new char[tempstr.length() + 1];
	char *bak = ptemp;
	if (ptemp == NULL)
	{
		printf("memory distribution failed...\n");
		return std::string::npos;
	}
	strcpy(ptemp, source.c_str());  //��string����ת��Ϊc��������
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
			continue; //һ��Ҫ�Ӳ�Ȼ���з��ʳ�ͻ
		}
	}
	*returnresult = bak; //��������
	//free(ptemp); //�����ٺ������ڲ��ͷţ���Ϊ����ͷžͻ���ʳ�ͻ,��main����free
	return 0; //����0����ɹ�
}



int find_terminal_codon(char *strcodon[], const char *gene_seq)  //������ֹ������
{
	if (strcodon == NULL)
	{
		return -1; //failure
	}
	for (char * *start = strcodon; start < strcodon + 3; start++)
	{
		std::string str1(gene_seq);
		size_t strlength = str1.length();
		int find_tail_string_pos = str1.rfind(*start); //��Ҫ�Ǻ�
		if (find_tail_string_pos == strlength - 3)
		{
			return 0;
		}
		else
		{
			continue; //û�ҵ��ͼ���
		}
	}
	return -1; //�����û�ҵ��ͷ���-1
}


std::string getGeneorgnism(std::string name_acquire) //��ȡfa�ļ���������
{
	std::string obtained_orgnism; //�洢��õ��ַ���
	if (name_acquire.empty()) //�����0����������,�����ݾ��Ǽٲ�ִ��������ز�����1����������
	{
		return "-1"; //failure ��gene_idΪ��
	}


	std::string match_str(":[[:blank:]]([[:alpha:]]+[[:blank:]][[:alpha:]]+)[[:blank:]]"); //������ʽ
	std::regex rule(match_str);
	std::smatch sm; //��ȡƥ����ַ���
	for (auto it = name_acquire.cbegin(); std::regex_search(it, name_acquire.cend(), sm, std::regex(match_str)); it = sm.suffix().first)
	{
		if (sm.size() != 0)
		{
			obtained_orgnism = sm.str();  //��ȡ����ð�ŵ�ƥ���ַ���
			break;
		}

	}

	return obtained_orgnism;
}


//��ʼ��ͷ�ڵ�
Queue *initialization()
{
	return NULL;
}

//ʵ�����ݵĲ���ṹ��
void enterQueuesort(Queue **phead, const char *id, std::string gene_seq)  //��Ӽ�����
{
	Queue *bak_for_return = NULL;
	//�����ǲ��ǿյ�ָ������������newһ���ڴ棬���ڴ��Ҫ����
	Queue *pnewelement = (Queue *)malloc(sizeof(Queue)); //��������
	//printf("queue=%p",pnewelement);
	pnewelement->gene_id = (char *)calloc(strlen(id) + 1, sizeof(char)); //�ֱ��ָ�뿪���ڴ�,���ҳ�ʼ��Ϊ0
	//printf("gene_id=%p\n", pnewelement->gene_id);
	pnewelement->gene_seq = (char *)calloc(gene_seq.length() + 1, sizeof(char));
	//�����ַ���
	strcpy(pnewelement->gene_id, id);
	strcpy(pnewelement->gene_seq, gene_seq.data());  //��return����
	pnewelement->pNext = NULL; //��һ��ָ�����NULL
	if (*phead == NULL) //���Ϊ�գ�ֱ�Ӽӵ��������(2019-4-23һֱ�����ԭ����������pheadû�и��Ǻţ�phead==NULL�Ǵ���ģ�*phead==NULL���ǶԵ�)
	{
		*phead = pnewelement;
		pnewelement->pNext = NULL;
	}
	else  //�����ͷ�ڵ��Ѿ���������
	{
		//�Ӵ�С������
		if (strlen(pnewelement->gene_seq) > strlen((*phead)->gene_seq))
		{
			pnewelement->pNext = *phead;
			*phead = pnewelement;
		}
		else
		{
			//β������
			Queue *ploop = *phead;
			while (ploop->pNext != NULL) //��ͷ��ѭ����β��
			{
				ploop = ploop->pNext;
			}
			if (strlen(pnewelement->gene_seq) <= strlen(ploop->gene_seq))
			{
				ploop->pNext = pnewelement;
				pnewelement->pNext = NULL; //ʵ��β����С�Ի�
			}
			//ʵ���в�����
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
				p1 = p1->pNext;//��ǰѭ��
			}
		}
	}
}


//�����������ڴ�
void enterqueue(Queue **phead, const char *gene_id, std::string gene_seq)
{

	if (*phead == NULL)
	{
		//Queue *pnewnode = (Queue *)calloc(1, sizeof(Queue));
		//Queue*pnewnode = (Queue *)malloc(sizeof(Queue));
		Queue *pnewnode = new Queue;
		if (pnewnode == NULL)
		{
			printf("�ڴ����ʧ��\n");
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
		//β��������
		Queue *ploop = *phead; //save head node
		while (ploop->pNext != NULL) //
		{
			ploop = ploop->pNext; //һֱѭ��
		}
		//Queue *pnewnode = (Queue *)malloc(sizeof(Queue));
		Queue *pnewnode1 = new Queue;
		if (pnewnode1 == NULL)
		{
			printf("�ڴ����ʧ��\n");
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
		ploop->pNext = pnewnode1; //������
	}
}

//�ͷ��ڴ�
Queue * freequeue(Queue *phead)
{
	Queue *p1, *p2; //��ֹҰָ��
	p1 = p2 = NULL;
	p1 = phead; //����ͷ�ڵ�

	//printf("p1->pNext=%s\n",p1->pNext->gene_id);
	//printf("p1->pNext->pNext=%s", p1->pNext->gene_seq);
	while (p1->pNext != NULL) //����Ϊʲôp1==NULL�ˣ���ôp1����˿�ָ���ˣ� //p1�ڵ㲻����p1�ṹ���pnext�����������Ľڵ㣬Ȼ��ɾ���м��
	{
		p2 = p1->pNext; //������һ���ڵ��ַ
		p1->pNext = p2->pNext; //�ͷŵ�p2�ڵ�   //������p1->pNext=p2->pNext����p1=p2->pNext,���˺þã�Ҫע����
		free(p2->gene_id);
		free(p2->gene_seq);
		delete p2; //�ͷŵ�p2�ڴ� 
		//printf("\n\n\n");
		//printQueue(phead); //�鿴ɾ�����м�״̬

	}
	free(phead->gene_id);
	free(phead->gene_seq);
	free(phead); //�ͷ�ͷ�ڵ�
	return NULL; //һ��Ҫ�����NULL
}


//��ʾ����
void printQueue(Queue *phead)
{
	using namespace std;
	const char * path = "C:\\Users\\xiaohui\\Desktop\\debug333.list";
	FILE *pout = fopen(path, "w");
	if (phead == NULL || pout == NULL)
	{
		printf("���ӡһ��������...\n");
		return;
	}
	while (phead != NULL)
	{
		fprintf(pout, "GeneID:%s\nGene_seq:%s\n", phead->gene_id, phead->gene_seq);
		cout << "geneID:" << phead->gene_id << "\n" << "GeneSeq:" << phead->gene_seq << "  length=" << strlen(phead->gene_seq) << endl;
		//pout << "geneID:" << phead->gene_id << "\n" << "GeneSeq:" << phead->gene_seq << "  length=" << strlen(phead->gene_seq) << "\n";
		//printf("%p,%p\n",phead,phead->pNext); //��ӡ�ڵ��ַ
		phead = phead->pNext;
	}
	fclose(pout);
	//return 1; //ִ�гɹ�����1
}

Queue* bubblesort(Queue *phead, char ch)  //ð�����򣬸��ݻ��򳤶�����
{
	if (phead == NULL || phead->pNext == NULL)
	{
		printf("your are in\n");
		printf("linktable empty or just a element!\n");
		return NULL;
	}
	else
	{
		if (ch == '>')   //�Ӵ�С����
		{
			for (Queue *p1 = phead; p1 != NULL; p1 = p1->pNext) //����ͷ�ڵ㣬ʹ�ñ���ѭ��
			{
				//int p1_len = strlen(p1->gene_seq); //���������õ����Ⱥ�Ƚϣ��Ǵ���ģ�������if����н��г����жϣ���������������ʾ
				//printf("p1_len=%d\n", p1_len);
				//printf("p1_seq=%s\n",p1->gene_seq);
				for (Queue *p2 = phead; p2 != NULL; p2 = p2->pNext)
				{

					//int p2_len = strlen(p2->gene_seq);
					//printf("p2_len=%d\n", p2_len);
					if (strlen(p1->gene_seq) > strlen(p2->gene_seq)) //�Ӵ�С����
					{ //ֻ��Ҫ����������ݴ洢���У������ǹ����
						Queue temp;
						temp.gene_id = p1->gene_id; //�����ָ�벻�ÿ������ַ������ÿ���
						temp.gene_seq = p1->gene_seq;

						p1->gene_id = p2->gene_id;
						p1->gene_seq = p2->gene_seq;
						p2->gene_id = temp.gene_id;
						p2->gene_seq = temp.gene_seq;
						//strcpy(p1->value, p2->value); //����ģ�������ʾ
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
			for (Queue *p1 = phead; p1 != NULL; p1 = p1->pNext) //����ͷ�ڵ㣬ʹ�ñ���ѭ��
			{
				//int p1_len = strlen(p1->gene_seq); //���������õ����Ⱥ�Ƚϣ��Ǵ���ģ�������if����н��г����жϣ���������������ʾ
				//printf("p1_len=%d\n", p1_len);
				//printf("p1_seq=%s\n",p1->gene_seq);
				for (Queue *p2 = phead; p2 != NULL; p2 = p2->pNext)
				{

					//int p2_len = strlen(p2->gene_seq);
					//printf("p2_len=%d\n", p2_len);
					if (strlen(p1->gene_seq) < strlen(p2->gene_seq)) //��С��������
					{ //ֻ��Ҫ����������ݴ洢���У������ǹ����
						Queue temp;
						temp.gene_id = p1->gene_id; //�����ָ�벻�ÿ������ַ������ÿ���
						temp.gene_seq = p1->gene_seq;

						p1->gene_id = p2->gene_id;
						p1->gene_seq = p2->gene_seq;
						p2->gene_id = temp.gene_id;
						p2->gene_seq = temp.gene_seq;
						//strcpy(p1->value, p2->value); //����ģ�������ʾ
						//strcpy(p1->key, p2->key);
						//strcpy(p2->value, temp.value);
						//strcpy(p2->key, temp.key);
					}
				}
			}
			return phead;
		}

	}
	//ʹ��ð������
	return NULL;
}


//��ȡ�ڵ�����Ҳ���ǻ��������
//дһ��ʵ���������¼�������������ʵ���ҵ����壬�����Լ��Լ��ķ����¼
int getnodenum(Queue *phead)
{
	int count = 0; //�ڵ������
	if (phead == NULL)
	{
		printf("������...\n");
		return 0;
	}
	Queue *pbak = phead;
	while (pbak != NULL)
	{
		count++;
		pbak = pbak->pNext;
	}
	return count;  //ͳ�ƽڵ����
}


//���������ض�Ƭ�β��ң��������ĳ��motif������ȫ���ҳ�����ʵ�֣� //��ʽ������������Ƕ��
void searchGeneMotif(Queue *phead, const char *motif, std::map<std::string, std::vector<std::string>> *myresult) //myresult�Ǵ�������
{  //������Ƕ��
	using namespace std;
	//multimap<string, string> *myresult = new multimap<string, string>;
	//multimap<string, string> myresult;
	string querystring(motif); //�����������motifת��Ϊstring����
	transform(querystring.begin(), querystring.end(), querystring.begin(), ::toupper); //ת��Ϊ��д
	if (phead == NULL)
	{
		printf("������...\n");
	}
	while (phead != NULL)
	{
		vector<std::string> vectemp;
		vector<std::string> tempvector;
		string subject(phead->gene_seq); //��char *ָ��ָ����ַ���ת��Ϊstring
		transform(subject.begin(), subject.end(), subject.begin(), ::toupper); //ת��Ϊ��д
		int count = 0;  //������
		int begin = -1; //��¼λ����Ϣ
		while ((begin = subject.find(querystring, begin + 1)) != string::npos)  //find�����Ƿ����ҵ��ĵ�һ����λ�ã���beginλ�ÿ�ʼ���ң�������0��
		{

			//insert(map<string,string>::value_type(cat_gene_id.data(),current->gene_seq));
			count++;
			char str[256] = { 0 };
			sprintf(str, "%d", begin + 1); //��1��ʼ��ӡ���û���
			string temp(str); //������ת��Ϊ�ַ���
			tempvector.push_back(temp);
			//vectemp.push_back(temp); //��λ����Ϣ�洢��vect������
			begin = begin + querystring.length();  //�ƹ����ҵ���һ���ַ����ĺ��棬���²���
		}
		if (count != 0) //˵����Gene seq�ҵ�motif
		{
			char str1[256] = { 0 };
			sprintf(str1, "%d", count);
			string temp1(str1);
			string geneSeq(phead->gene_seq);
			temp1 = "Number:" + temp1;
			//��λ����Ϣ��-->��ͷ���ӵ�һ�����һ���ַ���
			string catstring;
			for (size_t i = 1; i < tempvector.size(); ++i)
			{
				tempvector[0] += "-->";
				tempvector[0] += tempvector[i];
			}
			tempvector.resize(1); //���·����ڴ��С
			vectemp.push_back(tempvector[0]);
			vectemp.push_back(temp1); //��λ����Ϣ�洢��vect������
			vectemp.push_back(geneSeq);
			//outputsame_seqAndID.insert(map<string,string>::value_type(cat_gene_id.data(),current->gene_seq))
			//(*myresult).insert(make_pair(phead->gene_id, vectemp)); //�洢��������
			(*myresult).insert(map<string, vector<string>>::value_type(phead->gene_id, vectemp));
			//(*myresult).insert(pair<string, string>(phead->gene_id, phead->gene_seq));
			//(*myresult).insert(pair<string, string>(phead->gene_id, temp1)); //�������е�λ�ÿ�����ʱ��˳��
		}
		phead = phead->pNext; //ѭ��
	}
}


//���ҽڵ㣨���һ���
Queue *searchGeneID(Queue *phead, const char *key) //���ݻ�������
{
	if (phead == NULL)
	{
		printf("������...\n");
		return NULL;
	}
	char keyarray[1024] = { 0 };
	strcpy(keyarray, key);
	//keyarray[strlen(key) - 1] = '\0'; //ȥ��\n,�����Ǵ���ģ���Ѵ�����ַ�������������һ���ַ�
	for (size_t i = 0; i < strlen(keyarray); i++)  //Сдת��д
	{
		keyarray[i] = toupper(keyarray[i]); //26��Ӣ����ĸ��ת���������ַ�����
	}
	printf("searchkey=%s\n", keyarray);
	while (phead != NULL)
	{
		char pheadarray[1024] = { 0 };
		strcpy(pheadarray, phead->gene_id);
		printf("pheadarray=%s\n", pheadarray);
		for (size_t i = 0; i < strlen(pheadarray); i++)  //Сдת��д,��Ϊstrlen(p1)����Unsigned int,���㶨���iҲ������unsigned���߱����棬����ʹ��size_t
		{
			pheadarray[i] = toupper(pheadarray[i]); //26��Ӣ����ĸ��ת���������ַ�����
		}
		if (strcmp(keyarray, pheadarray) == 0) //�����ַ������
		{
			return phead; //���ص�ǰָ��ڵ�ĵ�ַ
		}
		phead = phead->pNext; //ѭ������
		//delete[]p1;  //�ҵ������˾Ͳ����ͷŶ��ڴ���
	}
	return NULL; //����������˻�û���ҵ�����NULL
}


//�޸�����
void modify(Queue *phead, const char *key, const char *value)
{
	if (phead == NULL)
	{
		return;
	}
	Queue *psearch = searchGeneID(phead, key);
	if (psearch == NULL)
	{
		printf("û�ҵ�..\n");
	}
	else
	{
		strcpy(psearch->gene_id, key);
		printf("�޸ĳɹ�..\n"); //���������û�����й�
	}
}


//�ַ�����ӣ��ַ���ɾ�����ַ�����
//�ж����������Ƿ���ȣ����߻���id�Ƿ��ظ���ʹ��strcmp��string.h��C���Ե�ͷ�ļ�

void deleteDuplicates(Queue *phead)
{
	using namespace std;
	//mapԪ��Ĭ�ϰ��ռ�����������
	map<string, string> outputsame_seqAndID; //�����洢����ַ����ļ���ֵ��gene_id��gene_seq�������Ǽ���ֵͬ��ͬ
	int count = 0; //��¼ɾ�����ظ�����
	int count1 = 0; //��¼id��ͬ������ͬ�ĸ���
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
				break; //ѭ������β������ѭ��
			}
			if (strcmp(current->gene_id, nextnode->gene_id) != 0 && nextnode && strcmp(current->gene_seq, nextnode->gene_seq) == 0) //�������һ������gene ID��һ����
			{
				string cat_gene_id(current->gene_id);
				cat_gene_id += "------>";
				cat_gene_id += nextnode->gene_id;
				outputsame_seqAndID.insert(map<string, string>::value_type(cat_gene_id.data(), current->gene_seq));
				count1++;
			}
			if (strlen(current->gene_seq) == strlen(nextnode->gene_seq) && strcmp(current->gene_id, nextnode->gene_id) == 0 && nextnode && strcmp(current->gene_seq, nextnode->gene_seq) == 0) //strcpm�������ַ�������0
			{
				count++;
				current->pNext = nextnode->pNext;
				//printf("current_seq=%s,nextnode=%s\n", current->gene_seq, nextnode->gene_seq);
				free(nextnode);
			}
			else
			{
				current = nextnode;  //ѭ��
			}
		}

	}
	if (!outputsame_seqAndID.empty()) //�����˲Ŵ����ļ�����ֹ���������ļ�,empty�������ǿշ���0���շ���1������ȡ��
	{
		const char *path1 = "C:\\Users\\xiaohui\\Desktop\\test_equal.fa"; //���ļ�����������Ŀ¼
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
	//�ۺ�����ǰ����
	fstream pf; //�������࣬ifstream,ofstream��fstream
	pf.open(path, ios::in); //��������ķ�ʽ��ȡ
	//�����
	fstream pfo;
	const char *outfile = "C:\\Users\\xiaohui\\Desktop\\outputfile2.txt";
	const char *outrevcom = "C:\\Users\\xiaohui\\Desktop\\revcomp.fa"; //������򻥲�����
	ofstream out(outrevcom, ios::out);
	//ofstream out(outfile, ios::out | ios::app); //��ͨд��򿪻���׷�Ӵ�//��ȷ��д��
	pfo.open(outfile, ios::out); //���ڴ�����д������ļ� //Ҳ��ȷ
	if (!pf.is_open() || !pfo.is_open() || !out.is_open())  //����Ǽٴ��ļ�ʧ��
	{
		cout << "Open file failure!\n";
		return -1;
	}

	//��ȡfa�ļ�����������Ŀǰ�޶���NCBI��fa�ļ���
	string getname;
	getline(pf, getname, '>'); //��ȡ��һ�е�gene_id,Ϊ��,����
	getline(pf, getname, '>'); //��ȡ�ڶ���
	size_t position_n = getname.find('\n'); //���ҵ�һ�λ��з����ֵ�λ�ã�Ϊ�˵õ������ID;size_t��unsigned int����
	getname.assign(getname, 0, position_n); //���ǶԵ�
	string orgnism_name_obtained = getGeneorgnism(getname);
	pf.clear(ios::goodbit); //�ļ�ָ������
	pf.seekg(ios::beg);
	//�����ʽ����,����
	pfo << "#Statistic File\t" << "Orgnism" << orgnism_name_obtained << endl; //ƥ����ַ����еõ�ð�źͿո�����Ͳ���Ҫð�źͿո�
	int gene_num = 0;
	//ͳ�ƻ������
	while (!pf.eof()) //û���ļ�ĩβ�ͼ���
	{
		char c;
		pf >> c;
		if (c == '>')
		{
			gene_num++;
		}
	}
	pf.clear(ios::goodbit);//��seekg֮ǰ��Ҫ�����������clear�����������ı�������������Ժ�Ϳ�����������seekg�����ˡ�
	pf.seekg(ios::beg); //begin����д,���ļ�ָ�����»ص�ͷ�� //seekp����������ļ�ָ��λ�ã�seekg�����������ļ�ָ��λ��
	/*ios::beg���ļ�������ʼλ��
	ios::cur���ļ����ĵ�ǰλ��
	ios::end���ļ����Ľ���λ��*/
	//����ļ���ʽ
	pfo << "#Gene_number: " << gene_num << endl;
	pfo << "#Gene_ID\t" << "Gene_length(bp)\t" << "A(bp)" << "\t" << "T(bp)" << "\t" << "G(bp)" << "\t" << "C(bp)" << "\t" << "N(bp)\t" << "GC(%)" << endl;


	string line;
	while (getline(pf, line, '>')) //����>�ָ����Ŷ�ȡ�ļ�
	{
		string gene_id; //��ȡ����id
		string gene_seq; //��ȡ��������
		int A, T, G, C, N; //ͳ�Ƽ����
		A = T = G = C = N = 0;
		double GC_content; //ͳ��GC����
		GC_content = 0.0;
		//if(line.find('\n')==string::npos);//npos����-1
		if (line.length() == 0 || line.find(0x0A) == 0)   //����������һ���Ƿ��ǿջ��߻س�,���з�����0λ�ã�0x0A=��\n��
		{
			continue; //������ͷ�Ŀ���
		}
		size_t pos_n = line.find('\n'); //���ҵ�һ�λ��з����ֵ�λ�ã�Ϊ�˵õ������ID;size_t��unsigned int����
		//gene_id = line.assign(line, 0, pos_n); �����
		gene_id.assign(line, 0, pos_n); //���ǶԵ�
		if (size_t pos_space = gene_id.find(0x20))  //�ҵ�space �ո�
		{
			gene_id.assign(gene_id, 0, pos_space); //��ȡ
		}
		if (size_t pos_tab = gene_id.find(0x09)) //�ҵ�ˮƽ�Ʊ��
		{
			gene_id.assign(gene_id, 0, pos_tab);
		}
		//gene_seq = line.substr(pos_n+1); //��ȡ���˻������Ƶ�ʣ������  //���������
		gene_seq.assign(line, pos_n + 1, line.length() - pos_n); //+1��Ϊ��ȥ����β�Ļ��з���  //���Ҳ������

		//remove������ͷ�ļ�#include<algorithm>�У����������ͷ�ļ�remove�����Ͳ���ʹ��
		//gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), '\n'), gene_seq.end()); 
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x0A), gene_seq.end()); //ɾ�����з���ͬ��
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x09), gene_seq.end()); //ɾ����ˮƽ�Ʊ����
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x08), gene_seq.end()); //ɾ��backspace
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x20), gene_seq.end()); //ɾ���ո�
		int length_seq = gene_seq.length(); //�������еĳ���,�������Ƿ�Ϊ3�ı���

		//�鿴���������Ƿ�ATG��ͷ  //transform��algorithmͷ�ļ��еĿ�
		transform(gene_seq.begin(), gene_seq.end(), gene_seq.begin(), ::toupper); //���ַ���ȫ��ת��Ϊ��д ��::tolower��Сд
		int atg_head = gene_seq.find("ATG"); //��ʼ������ //size_t��unsigned int
		//��ֹ������UAA��UGA��UAG
		char *termination_codon[3] = { "TAA","TGA","TAG" }; //ָ�����飬ÿһ��Ԫ��ָ��ָ��һ���ַ���(ͷָ����һ��������ַ�����һ��Ԫ��)

		int is_find_tail = find_terminal_codon(termination_codon, gene_seq.data());
		//���򻥲�����
		char *revp = NULL;
		if (reversecomplement(gene_seq, &revp) != string::npos) //����-1����ʧ��
		{
			string revcom(revp); //����������C�������ת��ΪC++�������
			string rev_gene = gene_id + "_ReverseComplementSeq";
			out << ">" << rev_gene << "\n" << revcom << "\n";  //��ӡ������ļ�
			free(revp); //�ͷŵ��ú���revcom���Ķ��ڴ�
		}

		//�ڴ�����֮��Ӽ�����
		//enterQueuesort(&phead, gene_id.data(), gene_seq); //���շ������Ķ��ڴ��׵�ַ(������)
		enterqueue(phead, gene_id.data(), gene_seq);
		//��ʾ����
		//printQueue(phead);
		//ͳ���ַ�
		for (string::iterator begin = gene_seq.begin(); begin != gene_seq.end(); begin++) //ʹ�õ�����ͳ���ַ���
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
		string::iterator end_find = gene_seq.end();  //���������������ڲ����ַ���

		//char *pbuf = (char *)calloc(gene_seq.length(),sizeof(char));//����Ҫ�����ﲻҪɾ��
		//strcpy(buf, gene_seq.data());  //data������ȡstring��c����ַ������� ����������ת����strcpy(buf,gene_seq.c_str);c_str��c�ַ���string�ļ�д
		//strcpy(pbuf, gene_seq.c_str());
		//free(pbuf); //������ͷ��ڴ�̫���������
		//fprintf(stderr,"c_str=%s\n",buf); //��ӡ�����������
		char gc_str[256] = { 0 };
		GC_content = ((G + C)*1.0 / (A + T + G + C)*1.0)*100.0;
		sprintf(gc_str, "%0.2f", GC_content);

		if (atg_head == 0 && length_seq % 3 == 0) //���ͷ����ATG���������ܱ�3�����ľͼ�һ���Ǻ�
		{
			gene_id += '*'; //��β�����һ���ַ�*
		}
		if (is_find_tail == 0) //���β��������ֹ�����ӵĻ������������Ǳ�ע
		{
			gene_id += '*';
		}
		pfo << gene_id << "\t" << length_seq << "\t" << A << "\t" << T << "\t" << G << "\t" << C << "\t" << N << "\t" << gc_str << endl;

	}
	pf.close();
	pfo.close();
	out.close(); //�ر��ļ���
	return 1;
}