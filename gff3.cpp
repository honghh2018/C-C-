#define _CRT_SECURE_NO_WARNINGS
#include "function.h"


 inline std::vector<std::string> split(std::string str, std::string pattern)  //字符串切割函数，返回一个vector数组（内联函数，加强函数的执行效率）
{
		std::string::size_type pos;
		std::vector<std::string> result;
		str += pattern;//扩展字符串以方便操作
		size_t size = str.size();
		for (size_t i = 0; i < size; i++)
		{
			pos = str.find(pattern, i);
			if (pos < size)
			{
				std::string s = str.substr(i, pos - i);
				result.push_back(s);
				i = pos + pattern.size() - 1;
			}
		}
	//NC_000913.3	RefSeq	gene	190	255	.	+	.	ID=gene-b0001;Dbxref=ASAP:ABE-0000006,ECOCYC:EG11277,EcoGene:EG11277,GeneID:944742;Name=thrL;gbkey=Gene;gene=thrL;gene_biotype=protein_coding;gene_synonym=ECK0001;locus_tag=b0001
	return result; //返回一个vector数组
}


 //get N number
 



void gff3statistic(const char *pathin,const char *pathout)
{
	time_t start, end;
	time(&start); //获取开始的时间
	using std::endl;
	unsigned long count_del_slash_n = 0;
	unsigned long gene_count = 0; //统计基因个数
	unsigned long transcript_count = 0; //统计转录本数
	unsigned long exon_count = 0; //统计外显子数
	unsigned long cds_count = 0;//统计cds
	std::ifstream in(pathin,std::ios::in);
	std::ofstream out(pathout, std::ios::out);
	const char *pathlog = "C:\\Users\\xiaohui\\Desktop\\gff3.statistic";
	std::ofstream outlog(pathlog, std::ios::out);
	if (!in.is_open() || !out.is_open())
	{
		printf("Your inputs files failed...\n");
		return;
	}
	outlog << "#gff3 statistic" << "\n";
	std::string line0;
	std::map<std::string,int> chrom;
	std::map<std::string, unsigned int> transnum;
	while (getline(in, line0))
	{
		if (line0.find("#") == 0) //跳过#开头的行
		{
			continue;
		}
		std::vector<std::string> mystr;
		mystr = split(line0, "\t");
		if (mystr[2] == "gene" || mystr[2] == "Gene" || mystr[2] == "GENE")
		{
			chrom[mystr[0]]++; //统计染色体数目和一条染色体上有多少个基因
			gene_count++;
		}
		if (mystr[2] == "mRNA" || mystr[2] == "mrna" || mystr[2] == "MRNA")
		{
			transnum[mystr[0]]++; //统计一个染色体上多少个转录本
			transcript_count++;
		}
		if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
		{
			exon_count++;
		}
		if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
		{
			cds_count++;
		}
	}
	outlog << "Gene numbers:" << gene_count << "\n"<<"Transcript numbers:"<< transcript_count<<"\n"<<"Exon numbers:"<<exon_count<<"\n"<<"CDS numbers:"<<cds_count<<endl;
	outlog << "Chromosome numbers:" << chrom.size() << endl; //打印染色体数目
	std::map<std::string, int>::iterator ibegin = chrom.begin();
	std::map<std::string, unsigned int>::iterator itranbegin = transnum.begin();
	for (; ibegin != chrom.end(); ++ibegin)
	{
			outlog << "Chromosome Name: " << ibegin->first << "\tGene numbers: " << ibegin->second << "\t"<<transnum[ibegin->first]<<endl;// "\t" << itranbegin->second << endl; //打印一条染色体上有多少条基因和多少个转录本
		//outlog << "Chromosome Name: " << ibegin->first << "\tGene numbers: " << ibegin->second << endl;// "\t" << itranbegin->second << endl; //打印一条染色体上有多少条基因和多少个转录本
			
	}
	in.clear(std::ios::goodbit);//在seekg之前需要调用流对象的clear方法，把流的标记清除掉，清除以后就可以正常调用seekg方法了。
	in.seekg(std::ios::beg); //begin的缩写,将文件指针重新回到头部 //seekp设置输出流文件指针位置，seekg设置输入流文件指针位置
	std::string line;
	while (getline(in, line)) //默认换行符为一行
	{
		if (line.find("#") == 0) //跳过#开头的行
		{
			continue;
		}
		std::vector<std::string> mystr;
		mystr = split(line, "\t");
		if (mystr[2] == "gene" || mystr[2] == "Gene" || mystr[2] == "GENE")
		{
			count_del_slash_n++;
			std::vector<std::string> mystr1= split(mystr[8], ";");
			std::vector<std::string> myname = split(mystr1[2], "=");
			out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << mystr1[0] << ":" << myname[1]<< endl;
		}
		if (mystr[2] == "mRNA" || mystr[2] == "mrna" || mystr[2] == "MRNA")
		{
			count_del_slash_n++;
			std::vector<std::string> mystr1 = split(mystr[8], ";");
			std::vector<std::string> myname = split(mystr1[3], "=");
			out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << mystr1[0] << ":" << myname[1] << ";" << mystr1[1]<< endl; //mystr1[1]是parent
		}
		if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
		{
			count_del_slash_n++;
			std::vector<std::string> mystr1 = split(mystr[8], ";");
			if (count_del_slash_n == (gene_count + transcript_count + exon_count + cds_count))
			{
				out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << mystr1[0];
			}
			else
			{
				out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << mystr1[0]<< ";" << mystr1[1]<< endl; //外显子没有Name
			}
			
		}
		if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
		{
			count_del_slash_n++;
			std::vector<std::string> mystr1 = split(mystr[8], ";");
			std::vector<std::string> myname = split(mystr1[3], "=");
			if (count_del_slash_n == (gene_count + transcript_count + exon_count + cds_count))
			{
				out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << mystr1[0] << ":" << myname[1]<<";" << mystr1[1];
			}
			else
			{
				out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << mystr1[0] << ":" << myname[1] << ";" << mystr1[1]<<endl;
			}
			
		}
//"NC_000913.3	RefSeq	gene	190	255	.	+	.	ID=gene-b0001;Dbxref=ASAP:ABE-0000006,ECOCYC:EG11277,EcoGene:EG11277,GeneID:944742;Name=thrL;gbkey=Gene;gene=thrL;gene_biotype=protein_coding;gene_synonym=ECK0001;locus_tag=b0001"
	}
	in.close();
	out.close();
	outlog.close();
	time(&end); //获取结束的时间
	printf("Elapse %d min\n", (unsigned int)(end - start) / 60); //打印时间
}

//正则表达式函数
inline bool is_regex(const std::string& object,std::vector<std::string> &myresult, const std::string &str)
{
	using std::cout;
	using std::endl;
	const std::regex pattern(str);
	std::smatch result;
	bool valid = std::regex_search(object, result, pattern);
	//此处result参数可有可无，result是一个字符串数组，用来存储正则表达式里面括号的内容。
	//std::sregex_iterator it(mrna.begin(), mrna.end(), pattern);
	//std::sregex_iterator end;
	//for (; it != end; ++it)
	//{
	//	cout << it->str() << endl;
	//}
	if (valid && result.size() > 0)
	{
		cout << result.size() << endl;
		for (size_t i = 1; i < result.size(); i++)
		{
			cout << result[i] << endl;
			myresult.push_back(result[i]);
		}
	}
		//cout << result.size() << endl;
		//for (size_t i = 1; i < result.size(); i++) //弃掉目标字符串
		//{
		//	cout << result[i]<<endl ;
		//	myresult.push_back(result[i]); //第一个是目标字符串本身，第二个是第一个小括号，依此类推
	//	}
	
	return valid;
}

//string 字符串转换为数字后排序
std::vector<int> sortstring(std::vector<std::string> &in)
{
	std::vector<int> myreturn;
	printf("i am in\n");
	for (size_t i = 0; i < in.size(); i++)
	{
		myreturn.push_back(atoi(in[i].c_str()));
	}
	std::sort(myreturn.begin(), myreturn.end(), std::less<int>()); //升序排序
	return myreturn;
}

inline std::string tr_reverse(std::string &str)
{
	size_t len = str.length();
	char *pstr = new char[len];
	strcpy(pstr, str.c_str());
	for (size_t i = 0; i < len; i++)
	{
		if (pstr[i] == 'A')
		{
			pstr[i] = 'T';
		}
		else if (pstr[i] == 'T')
		{
			pstr[i] = 'A';
		}
		else if (pstr[i] == 'C')
		{
			pstr[i] = 'G';
		}
		else if (pstr[i] == 'G')
		{
			pstr[i] = 'C';
		}
		else if (pstr[i] == 'R')
		{
			pstr[i] = 'Y';
		}
		else if (pstr[i] == 'Y')
		{
			pstr[i] = 'R';
		}
		else if (pstr[i] == 'M')
		{
			pstr[i] = 'K';
		}
		else if (pstr[i] == 'K')
		{
			pstr[i] = 'M';
		}
		else if (pstr[i] == 'H')
		{
			pstr[i] = 'D';
		}
		else if (pstr[i] == 'D')
		{
			pstr[i] = 'H';
		}
		else if (pstr[i] == 'B')
		{
			pstr[i] = 'V';
		}
		else if (pstr[i] == 'V')
		{
			pstr[i] = 'B';
		}
	}
	//$trans_seq = ~tr / RYMKHDBV / YRKMDHVB / ;
		int tmp;
		for (size_t i = 0; i < len / 2; i++)
		{
			tmp = pstr[i];
			pstr[i] = pstr[len - 1 - i];
			pstr[len - 1 - i] = tmp;
		}
		std::string myre(pstr);
		return myre;
		delete[]pstr;
}

//check gff3 whether or not it had LncRNA, rRNA,tRNA,CircRNA or MiRNA
int check_no_coding_gene(const char *ingff3,size_t &linenum)   //linenum用于统计文件数量
{
	enum {correct=1,incorrect=0,rRNA1=2,tRNA1=3,lncRNA1=4,miRN1A=5,circRNA1=6,noexon=7};  //定义匿名枚举类型
	short markgene=0, markmrna=0, markexon=0, markcds=0,rRNA=0,tRNA=0,LncRNA=0,CircRNA=0,MiRNA=0;
	int rRNA_count = 0, tRNA_count = 0, circRNA_count = 0, LncRNA_count = 0, MiRNA_count = 0, cds_count = 0, gene_count = 0,mRNA_count=0;
	short no_coding_rna_exon_count= 0;
	std::map<std::string, size_t>stat_exon;
	std::map<std::string, size_t> stat_rna;
	std::map<std::string, size_t> stat_gene;
	//rRNA匹配
	//boost::regex rRNA_pattern("R[ibosomal\\s]{0,9}RNA", boost::regex::perl | boost::regex::icase);
	//boost::smatch rRNA_result; //不用报错捕获的结果，就不用定义这个变量
	boost::smatch markresult;
	std::string xrna;
	std::ifstream in(ingff3, std::ios::in);
	if (!in.is_open())
	{
		std::cout << "Open file *.gff3 failed." << std::endl;
		return incorrect; //返回一个枚举类型，代表的常量值是0，default为int类型
	}
	std::string line;
	while (getline(in, line))
	{
		if (line.find("#") == 0 || line.length() == 0 || line.find(0x0A) == 0) //跳过#开头的行 //用来查找上一行是否是空或者回车,换行符都在0位置，0x0A=‘\n’
		{
			continue;
		}
		std::vector<std::string> mystr(split(line,"\t")); //保存切割的字符串
		if (boost::regex_match(mystr[2], markresult, boost::regex("(gene)", boost::regex::perl | boost::regex::icase)))
		{
			gene_count++;
			markgene = 1;
			continue;

		}
		//xrna_result 要事先定义
		if (boost::regex_match(mystr[2], markresult,boost::regex("([\\w_-\\s]{0,15}RNA)", boost::regex::perl | boost::regex::icase))) //匹配到*RNA,包括rRNA,tRNA，LncRNA等等
		{
			markmrna = 1;
			stat_rna[markresult[1]]++; //统计不同的mRNA数量
		}
		if (boost::regex_match(mystr[2], markresult, boost::regex("(exon)", boost::regex::perl | boost::regex::icase))) //匹配exon
		{

		}
		if (boost::regex_match(mystr[2], markresult, boost::regex("(cds)", boost::regex::perl | boost::regex::icase)))
		{

		}
		/*
		if (mystr[2] == "mRNA" || mystr[2] == "mrna" || mystr[2] == "MRNA")
		{
			mRNA_count++;
			xrna = mystr[2];
			markmrna = 1;
			continue;

		}
		*/
		//match是完全匹配，search是搜索子串
		if (boost::regex_match(mystr[2],boost::regex("R[ibosomal_-\\s]{0,10}RNA", boost::regex::perl | boost::regex::icase))) //核糖体RNA，只有exon没有cds ,不要将匹配的结果放到一个变量里就不用第二个参数result
		{
			xrna = mystr[2];
			rRNA_count++;
			rRNA = 1;
			continue;
		}
		if (boost::regex_match(mystr[2], boost::regex("T[ransfer\\s_-ribonucleic_-\\sacid_-\\sRNA]{0,34}", boost::regex::perl | boost::regex::icase))) //转运RNA只有exon没有cds
		{
			xrna = mystr[2];
			tRNA_count++;
			tRNA = 1;
			continue;
		}
		if (boost::regex_match(mystr[2], boost::regex("l[ong_-\\snoncoding_-\\s]{0,19}RNA", boost::regex::perl | boost::regex::icase))) //长链非编码
		{
			xrna = mystr[2];
			LncRNA_count++;
			LncRNA = 1;
			continue;
		}
		if (boost::regex_match(mystr[2], boost::regex("Mi[cro_-\\s]{0,9}RNA", boost::regex::perl | boost::regex::icase))) //小RNA
		{
			xrna = mystr[2];
			MiRNA_count++;
			MiRNA = 1;
			continue;
		}
		if (boost::regex_match(mystr[2], boost::regex("c[irc_-\\s]{0,7}RNA", boost::regex::perl | boost::regex::icase)) ) //circRNA
		{
			xrna = mystr[2];
			circRNA_count++;
			CircRNA = 1;
			continue;
		}
		if ((mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON"))
		{
			stat_exon[xrna]++;  //统计不同rna的exon个数
			markexon = 1;
			continue;
		}
		if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
		{
			cds_count++;
			markcds = 1;
			continue;

		}
	}
	if (!(markexon==1)) //如果没有外显子，也就是只有cds
	{
		std::cout << "gff3 have no exon mark.\n" << std::endl;
		if (markgene &&markmrna &&markcds) //有gene和mrna和cds
		{
			linenum = cds_count + mRNA_count + gene_count;
			return noexon;
		}
		if (markgene &&markcds) //有gene和cds
		{
			linenum = cds_count +gene_count;
			return noexon;
		}
		if (markmrna &&markcds) //有mrna和cds
		{
			linenum = cds_count + mRNA_count;
			return noexon;
		}
	}

	//如果有外显子
	if (markgene && markmrna && markexon && markcds &&CircRNA &&MiRNA &&LncRNA &&rRNA &&tRNA) //这些全部都有
	{
		std::map<std::string, size_t>::iterator begin = stat_exon.begin();
		for (; begin != stat_exon.end(); begin++)
		{
			linenum += begin->second;
		}
		linenum += cds_count + circRNA_count + MiRNA_count + LncRNA_count + rRNA_count + tRNA_count + mRNA_count + gene_count; //如果没有的就是0
		return 1;
	}
	//
	if (markgene && markmrna && markexon && markcds &&CircRNA &&MiRNA &&rRNA&&tRNA) //只有crna,rrna,mirna,rRNA
	{
		std::map<std::string, size_t>::iterator begin = stat_exon.begin();
		for (; begin != stat_exon.end(); begin++)
		{
			linenum += begin->second;
		}
		linenum += cds_count + circRNA_count + MiRNA_count +rRNA_count + tRNA_count + mRNA_count + gene_count;
		return 1;
	}
	if (markgene && markmrna && markexon && markcds &&CircRNA &&MiRNA &&rRNA&&tRNA) //只有crna,rrna,mirna,rRNA
	{
		std::map<std::string, size_t>::iterator begin = stat_exon.begin();
		for (; begin != stat_exon.end(); begin++)
		{
			linenum += begin->second;
		}
		linenum += cds_count + circRNA_count + MiRNA_count + rRNA_count + tRNA_count + mRNA_count + gene_count;
		return 1;
	}
	if (markgene && markmrna && markexon && markcds &&CircRNA) //只有circrna
	{
		std::map<std::string, size_t>::iterator begin = stat_exon.begin();
		for (; begin != stat_exon.end(); begin++)
		{
			linenum += begin->second;
		}
		linenum += cds_count + circRNA_count + mRNA_count + gene_count; //如果没有的就是0
		return 1;
	}
	if (markgene && markmrna && markexon && markcds &&MiRNA) //只有MiRNA
	{
		std::map<std::string, size_t>::iterator begin = stat_exon.begin();
		for (; begin != stat_exon.end(); begin++)
		{
			linenum += begin->second;
		}
		linenum += cds_count + MiRNA_count + mRNA_count + gene_count; //如果没有的就是0
		return 1;
	}
	if (markgene && markmrna && markexon && markcds &&rRNA) //只有rrna
	{
		std::map<std::string, size_t>::iterator begin = stat_exon.begin();
		for (; begin != stat_exon.end(); begin++)
		{
			linenum += begin->second;
		}
		linenum += cds_count + rRNA_count + mRNA_count + gene_count; //如果没有的就是0
		return 1;
	}
	if (markgene && markmrna && markexon && markcds && LncRNA) //只有lncrna
	{
		std::map<std::string, size_t>::iterator begin = stat_exon.begin();
		for (; begin != stat_exon.end(); begin++)
		{
			linenum += begin->second;
		}
		linenum += cds_count + LncRNA_count + mRNA_count + gene_count; //如果没有的就是0
		return 1;
	}
	if (markgene && markmrna && markexon && markcds && tRNA) //只有trna
	{
		std::map<std::string, size_t>::iterator begin = stat_exon.begin();
		for (; begin != stat_exon.end(); begin++)
		{
			linenum += begin->second;
		}
		linenum += cds_count + tRNA_count + mRNA_count + gene_count; //如果没有的就是0
		return 1;
	}
	//没有gene
	if (markmrna && markexon && markcds &&CircRNA &&MiRNA &&LncRNA &&rRNA &&tRNA)
	{

	}
	if (markmrna && markexon && markcds && tRNA) //转运rna
	{

	}
}









int checkgff3mark(const char *pathingff3, const char *outgff3)  //生成一个符合下面正则表达式的gff3文件作为下面函数的输入
{
	int markgene = 0, flaggene = 0;
	int markmrna = 0, flagmrna = 0;
	int markexon = 0;
	int markcds = 0;
	int count = 1000; //读取1000行判断文件是否是标准的gff3文件
	std::ifstream in(pathingff3, std::ios::in);
	std::string line;
	boost::regex pattern("ID=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)", boost::regex::perl | boost::regex::icase); //++是流模式减少回溯，提高正则速度
	boost::smatch result;
	boost::regex pattern1("ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)[^\\s\"]+gene=([\\w\\d]++);", boost::regex::perl | boost::regex::icase); //++是流模式减少回溯，提高正则速度
	boost::smatch result1;
	boost::regex pattern2("ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)", boost::regex::perl | boost::regex::icase);
	boost::smatch result2;
	//mrna 的匹配
	boost::regex pattern_mrna("ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)", boost::regex::perl | boost::regex::icase);
	boost::smatch result_mrna;
	boost::regex pattern_mrna1("ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+", boost::regex::perl | boost::regex::icase); //没有name的匹配
	boost::smatch result_mrna1;
	//gene only id
	boost::regex geneid("ID=([^;\\s\"]++)[^\"]+", boost::regex::perl | boost::regex::icase);
	boost::smatch gene_id;
	//exon 匹配
	boost::regex exon_id("ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]++)", boost::regex::perl | boost::regex::icase);
	boost::smatch exon_result;
	//匹配cds
	boost::regex cds_match("ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)", boost::regex::perl | boost::regex::icase);
	boost::smatch cds_match_result;
	boost::regex cds_match1("ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+", boost::regex::perl | boost::regex::icase);
	boost::smatch cds_match_result1;
	//没有ID的外显子,Name就是它的ID
	boost::regex exon_no_id("Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)", boost::regex::perl | boost::regex::icase);
	boost::smatch exon_no_id_result;
	std::map<std::string, std::string> gene;
	std::map<std::string, std::map<std::string, std::string>> mrnaname;
	std::map< std::string, std::map<std::string, std::string >> cdsname;
	std::cout << "markgene=" << markgene << "\tmarkmrna=" << markmrna << "\tmarkexon=" << markexon << "\tmarkcds=" << markcds << std::endl;

	while (getline(in, line))
	{

		if (line.find("#") == 0 || line.length() == 0 || line.find(0x0A) == 0) //跳过#开头的行 //用来查找上一行是否是空或者回车,换行符都在0位置，0x0A=‘\n’
		{
			continue;
		}
		if (count == 0)
		{
			break; //跳出循环
		}
		count--; //循环50次//1000
		std::vector<std::string> mystr(split(line, "\t")); //保存切割的字符串
		if (mystr[2] == "gene" || mystr[2] == "Gene" || mystr[2] == "GENE")
		{
			markgene = 1;
			continue;

		}
		if (mystr[2] == "mRNA" || mystr[2] == "mrna" || mystr[2] == "MRNA")
		{
			markmrna = 1;
			continue;

		}
		if (mystr[2] == "rRNA" || mystr[2] == "rrna" || mystr[2] == "RRNA") //核糖体RNA，只有exon没有cds
		{

		}
		if (mystr[2] == "tRNA" || mystr[2] == "trna" || mystr[2] == "TRNA") //转运RNA只有exon没有cds
		{

		}
		if (mystr[2] == "LncRNA" || mystr[2] == "lncrna" || mystr[2] == "LNCRNA") //长链非编码
		{

		}
		//if()
		if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
		{
			markexon = 1;
			continue;
		}
		if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
		{
			markcds = 1;
			continue;

		}


	}

	std::cout << "check gff3 completion..." << std::endl;
	std::cout << "markgene=" << markgene << "\tmarkmrna=" << markmrna << "\tmarkexon=" << markexon << "\tmarkcds=" << markcds << std::endl;
	in.clear(std::ios::goodbit);//在seekg之前需要调用流对象的clear方法，把流的标记清除掉，清除以后就可以正常调用seekg方法了。
	in.seekg(std::ios::beg); //begin的缩写,将文件指针重新回到头部 //seekp设置输出流文件指针位置，seekg设置输入流文件指针位置
	std::ofstream out(outgff3, std::ios::out);
	std::string line1;
	int flag1 = 0;
	std::string SaveGenenName;
	std::string SaveGeneID;
	while (getline(in, line1))
	{
		if (line1.find("#") == 0 || line1.length() == 0 || line1.find(0x0A) == 0) //跳过#开头的行 //用来查找上一行是否是空或者回车,换行符都在0位置，0x0A=‘\n’
		{
			continue;
		}
		std::vector<std::string> mystr(split(line1, "\t")); //保存切割的字符串
		if (markgene && markmrna && markexon && markcds)
		{
			//printf("i am in\n");
			if (mystr[2] == "gene" || mystr[2] == "Gene" || mystr[2] == "GENE")
			{
				//"ID=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)", 
				if (boost::regex_search(mystr[8], result, pattern))
				{
					//printf("i am in sec...\n");
					SaveGenenName = result[2];
					flag1 = 1;
					///std::cout << result[0] << std::endl;
					//std::cout << result[1] << std::endl;
					//std::cout << result[2] << std::endl;
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result[1] << ";" << "Name=" << result[2] << std::endl;
					continue;
				}
				//ID=([^;\\s\"]++)[^\"]+  只有id没有name的gene
				if (boost::regex_search(mystr[8], gene_id, geneid))
				{
					//printf("i am in third...\n");
					SaveGeneID = gene_id[1];
					flag1 = 2;
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << gene_id[1] << ";" << "Name=" << gene_id[1] << std::endl; //gene没有name的直接用id代替name
				}
				continue;
			}
			if (mystr[2] == "mRNA" || mystr[2] == "mrna" || mystr[2] == "MRNA")
			{
				//printf("mrna in...\n");
				//"ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)[^\\s\"]+gene=([\\w\\d]++);",
				if (boost::regex_search(mystr[8], result1, pattern1))
				{

					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result1[1] << ";" << "Parent=" << result[1] << ";Name=" << result1[2] << ";gene=" << result1[3] << std::endl;
					continue;
				}
				//"ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)" 没有Name
				if (boost::regex_search(mystr[8], result2, pattern2))
				{
					//std::cout <<"result2="<< result2[0] << std::endl;
					//std::cout << "result2="<<result2[1] << std::endl;
					//std::cout << "result2="<<result2[2] << std::endl;
					//std::cout << "result="<<result[2] << std::endl;
					if (flag1 == 1)
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result2[1] << ";Parent=" << result2[2] << ";Name=" << result2[1] << ";gene=" << SaveGenenName << std::endl; //有name就等于gene name
						flag1 = 0;
						continue;
					}
					if (flag1 == 2)
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result2[1] << ";Parent=" << result2[2] << ";Name=" << result2[1] << ";gene=" << SaveGeneID << std::endl; //没name就等于gene id
						flag1 = 0;
						continue;
					}
				}
				continue;
			}
			if (in.eof()) //判断如果是读取的最后一行，那么输入不打印换行
			{
				if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
				{

					//
					if (boost::regex_search(mystr[8], exon_result, exon_id))
					{
						printf("i am in exon...\n");
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << exon_result[1] << ";Parent=" << exon_result[2];
						continue;
					}
					//exon_no_id //"Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)", boost::regex::perl | boost::regex::icase;
					if (boost::regex_search(mystr[8], exon_no_id_result, exon_no_id))
					{
						printf("i am in exon111...\n");
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << exon_no_id_result[2] << ";Parent=" << exon_no_id_result[1];
						continue;
					}
					continue;
				}
				if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
				{
					if (boost::regex_search(mystr[8], cds_match_result, cds_match))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3]; //cds有name就用name
						continue;
					}
					if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1]; //没有Name的就等于cds的ID
						continue;
					}
				}
			} //如果是最后一行不打印换行符

			if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
			{

				//
				if (boost::regex_search(mystr[8], exon_result, exon_id))
				{
					printf("i am in exon...\n");
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << exon_result[1] << ";Parent=" << exon_result[2] << std::endl;
					continue;
				}
				//exon_no_id //"Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)", boost::regex::perl | boost::regex::icase;
				if (boost::regex_search(mystr[8], exon_no_id_result, exon_no_id))
				{
					printf("i am in exon111...\n");
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << exon_no_id_result[2] << ";Parent=" << exon_no_id_result[1] << std::endl;
					continue;
				}
				continue;
			}
			if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
			{
				if (boost::regex_search(mystr[8], cds_match_result, cds_match))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3] << std::endl; //cds有name就用name
					continue;
				}
				if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1] << std::endl; //没有Name的就等于cds的ID
					continue;
				}
			}
			continue; //如果执行这个语句就跳过下面的语句
		}
		if (markgene && markcds &&markexon) //有gene,cds和exon，没有mrna执行这个语句
		{
			if (mystr[2] == "gene" || mystr[2] == "Gene" || mystr[2] == "GENE")
			{
				//"ID=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)", 
				if (boost::regex_search(mystr[8], result, pattern))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result[1] << ";" << "Name=" << result[2] << std::endl;
					out << mystr[0] << "\t" << mystr[1] << "\t" << "mRNA" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=mrna:" << result[1] << "_rep" << ";" << "Parent=" << result[1] << ";Name=" << "mrna:" << result[1] << "_rep" << ";gene=" << result[2] << std::endl;
					continue;
				}
				//ID=([^;\\s\"]++)[^\"]+  只有id没有name的gene
				if (boost::regex_search(mystr[8], gene_id, geneid))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << gene_id[1] << ";" << "Name=" << gene_id[1] << std::endl; //gene没有name的直接用id代替name
					out << mystr[0] << "\t" << mystr[1] << "\t" << "mRNA" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=mrna:" << result[1] << "_rep" << ";" << "Parent=" << result[1] << ";Name=" << "mrna:" << result[1] << "_rep" << ";gene=" << result[1] << std::endl; //gene没有name的直接用id代替name
				}
				continue;
			}
			if (in.eof()) //判断如果是读取的最后一行，那么输入不打印换行
			{
				if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
				{
					if (boost::regex_search(mystr[8], cds_match_result, cds_match))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3]; //cds有name就用name
						continue;
					}
					if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1]; //没有Name的就等于cds的ID
						continue;
					}
				}
				if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
				{
					if (boost::regex_search(mystr[8], exon_result, exon_id))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << exon_result[1] << ";Parent=" << exon_result[2];
					}
					continue;
				}
				continue;
			} //如果是末尾了，不添加换行

			if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
			{
				if (boost::regex_search(mystr[8], cds_match_result, cds_match))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3] << std::endl; //cds有name就用name
					continue;
				}
				if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1] << std::endl; //没有Name的就等于cds的ID
					continue;
				}
			}
			if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
			{
				if (boost::regex_search(mystr[8], exon_result, exon_id))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << exon_result[1] << ";Parent=" << exon_result[2] << std::endl;
				}
				continue;
			}
			continue;
		}
		if (markgene &&markmrna&& markcds) //有gene,mrna和exon，没有cds执行这个
		{
			if (mystr[2] == "gene" || mystr[2] == "Gene" || mystr[2] == "GENE")
			{
				//"ID=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)", 
				if (boost::regex_search(mystr[8], result, pattern))
				{
					flag1 = 1;
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result[1] << ";" << "Name=" << result[2] << std::endl;
					continue;
				}
				//ID=([^;\\s\"]++)[^\"]+  只有id没有name的gene
				if (boost::regex_search(mystr[8], gene_id, geneid))
				{
					flag1 = 2;
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << gene_id[1] << ";" << "Name=" << gene_id[1] << std::endl; //gene没有name的直接用id代替name
				}
				continue;
			}
			if (mystr[2] == "mRNA" || mystr[2] == "mrna" || mystr[2] == "MRNA")
			{
				//"ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)[^\\s\"]+gene=([\\w\\d]++);",
				if (boost::regex_search(mystr[8], result1, pattern1))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result1[1] << ";" << "Parent=" << result[1] << ";Name=" << result1[2] << ";gene=" << result1[3] << std::endl;
					continue;
				}
				//"ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)" 没有Name
				if (boost::regex_search(mystr[8], result2, pattern2))
				{
					if (flag1 == 1)
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result2[1] << ";Parent=" << result2[2] << ";Name=" << result2[1] << ";gene =" << result[2] << std::endl; //有name就等于gene name
						flag1 = 0;
						continue;
					}
					if (flag1 == 2)
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result2[1] << ";Parent=" << result2[2] << ";Name=" << result2[1] << ";gene =" << result[1] << std::endl; //没name就等于gene id
						flag1 = 0;
						continue;
					}
				}
				continue;
			}
			if (in.eof()) //如果读到文件末尾了，不换行
			{
				if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
				{
					//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+
					if (boost::regex_search(mystr[8], exon_result, exon_id))
					{
						if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 0)
						{
							out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "0" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
						}
						if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 1)
						{
							out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "1" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
						}
						if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 2)
						{
							out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "2" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
						}
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << exon_result[1] << ";Parent=" << exon_result[2]; //不换行
						continue;
					}
					continue;
				}
			}//不换行

			if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
			{
				//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+
				if (boost::regex_search(mystr[8], exon_result, exon_id))
				{
					if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 0)
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "0" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
					}
					if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 1)
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "1" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
					}
					if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 2)
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "2" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
					}
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << exon_result[1] << ";Parent=" << exon_result[2] << std::endl;
					continue;
				}
				continue;
			}

			continue;
		}
		if (markmrna &&markcds &&markexon) //有mrna,exon和cds没有gene
		{
			if (mystr[2] == "mRNA" || mystr[2] == "mrna" || mystr[2] == "MRNA")
			{
				//([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)"; //这里没有parent，因为自己就是带头大哥
				//ID=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)
				if (boost::regex_search(mystr[8], result, pattern))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << "gene" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=gene:" << result[1] << "_rep" << ";Name=gene:" << result[1] << "_rep" << std::endl;
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result[1] << ";" << "Parent=gene:" << result[1] << "_rep" << ";Name=" << result[2] << ";gene=gene:" << result[1] << "_rep" << std::endl;
					continue;
				}
				//"ID=([^;\\s\"]++)[^\"]+" 没有Name
				if (boost::regex_search(mystr[8], gene_id, geneid))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << "gene" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=gene:" << gene_id[1] << "_rep" << ";Name=gene:" << gene_id[1] << "_rep" << std::endl; //name是自己的ID
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << gene_id[1] << ";Parent=gene:" << gene_id[1] << "_rep" << ";Name=" << gene_id[1] << ";gene =gene:" << gene_id[1] << "_rep" << std::endl; //name为自己的ID
					continue;
				}
				continue;
			}
			if (in.eof()) //如果读到文件末尾了，不换行
			{
				if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
				{
					if (boost::regex_search(mystr[8], cds_match_result, cds_match))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3]; //cds有name就用name
						continue;
					}
					if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1]; //没有Name的就等于cds的ID
						continue;
					}
				}
				if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
				{
					//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+
					if (boost::regex_search(mystr[8], exon_result, exon_id))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << exon_result[1] << ";Parent=" << exon_result[2];
					}
					continue;
				}
			} //不换行
			if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
			{
				if (boost::regex_search(mystr[8], cds_match_result, cds_match))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3] << std::endl; //cds有name就用name
					continue;
				}
				if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1] << std::endl; //没有Name的就等于cds的ID
					continue;
				}
			}
			if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
			{
				//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+
				if (boost::regex_search(mystr[8], exon_result, exon_id))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << exon_result[1] << ";Parent=" << exon_result[2] << std::endl;
				}
				continue;
			}

			continue;
		}
		//到这里
		if (markmrna &&markcds) //有mrna和cds没有gene和exon
		{
			if (mystr[2] == "mRNA" || mystr[2] == "mrna" || mystr[2] == "MRNA")
			{
				//([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)"; //这里没有parent，因为自己就是带头大哥
				//ID=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)
				if (boost::regex_search(mystr[8], result, pattern))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << "gene" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=gene:" << result[1] << "_rep" << ";Name=gene:" << result[1] << "_rep" << std::endl;
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result[1] << ";" << "Parent=gene:" << result[1] << "_rep" << ";Name=" << result[2] << ";gene=gene:" << result[1] << "_rep" << std::endl;
					continue;
				}
				//"ID=([^;\\s\"]++)[^\"]+" 没有Name
				if (boost::regex_search(mystr[8], gene_id, geneid))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << "gene" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=gene:" << gene_id[1] << "_rep" << ";Name=gene:" << gene_id[1] << "_rep" << std::endl; //name是自己的ID
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << gene_id[1] << ";Parent=gene:" << gene_id[1] << "_rep" << ";Name=" << gene_id[1] << ";gene =gene:" << gene_id[1] << "_rep" << std::endl; //name为自己的ID
					continue;
				}
				continue;
			}
			if (in.eof()) //如果读到文件末尾了，不换行
			{
				if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
				{
					//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)
					if (boost::regex_search(mystr[8], cds_match_result, cds_match))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3] << std::endl; //cds有name就用name
						out << mystr[0] << "\t" << mystr[1] << "\t" << "exon" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "." << "\tID=exon:" << cds_match_result[1] << "_rep" << ";Parent=" << cds_match_result[2]; //cds有name就用name

						continue;
					}
					//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+
					if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1] << std::endl; //没有Name的就等于cds的ID
						out << mystr[0] << "\t" << mystr[1] << "\t" << "exon" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "." << "\tID=exon:" << cds_match_result1[1] << "_rep" << ";Parent=" << cds_match_result1[2]; //没有Name的就等于cds的ID
						continue;
					}
				}
			} //末尾不换行

			if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
			{
				//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)
				if (boost::regex_search(mystr[8], cds_match_result, cds_match))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3] << std::endl; //cds有name就用name
					out << mystr[0] << "\t" << mystr[1] << "\t" << "exon" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "." << "\tID=exon:" << cds_match_result[1] << "_rep" << ";Parent=" << cds_match_result[2] << std::endl; //cds有name就用name

					continue;
				}
				//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+
				if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1] << std::endl; //没有Name的就等于cds的ID
					out << mystr[0] << "\t" << mystr[1] << "\t" << "exon" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "." << "\tID=exon:" << cds_match_result1[1] << "_rep" << ";Parent=" << cds_match_result1[2] << std::endl; //没有Name的就等于cds的ID
					continue;
				}
			}

			continue;
		}
		if (markgene &&markcds) //有gene和cds没有mrna和exon
		{
			if (mystr[2] == "gene" || mystr[2] == "Gene" || mystr[2] == "GENE")
			{
				//"ID=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)", 
				if (boost::regex_search(mystr[8], result, pattern))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result[1] << ";" << "Name=" << result[2] << std::endl;
					out << mystr[0] << "\t" << mystr[1] << "\t" << "mRNA" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=mrna:" << result[1] << "_rep" << ";" << "Parent=" << result[1] << ";Name=" << "mrna:" << result[1] << "_rep" << ";gene=" << result[2] << std::endl;
					continue;
				}
				//ID=([^;\\s\"]++)[^\"]+  只有id没有name的gene
				if (boost::regex_search(mystr[8], gene_id, geneid))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << gene_id[1] << ";" << "Name=" << gene_id[1] << std::endl; //gene没有name的直接用id代替name
					out << mystr[0] << "\t" << mystr[1] << "\t" << "mRNA" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=mrna:" << result[1] << "_rep" << ";" << "Parent=" << result[1] << ";Name=" << "mrna:" << result[1] << "_rep" << ";gene=" << result[1] << std::endl; //gene没有name的直接用id代替name
				}
				continue;

			}
			if (in.eof()) //末尾不换行
			{
				if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
				{
					//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)
					if (boost::regex_search(mystr[8], cds_match_result, cds_match))
					{
						//printf("i am in cds...\n");
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3] << std::endl; //cds有name就用name
						out << mystr[0] << "\t" << mystr[1] << "\t" << "exon" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "." << "\t" << "ID=exon:" << cds_match_result[1] << "_rep" << ";Parent=" << cds_match_result[2]; //cds有name就用name

						continue;
					}
					//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+
					if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
					{
						//printf("i am in cds1...\n");
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1] << std::endl; //没有Name的就等于cds的ID
						out << mystr[0] << "\t" << mystr[1] << "\t" << "exon" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "." << "\t" << "ID=exon:" << cds_match_result1[1] << "_rep" << ";Parent=" << cds_match_result1[2]; //没有Name的就等于cds的ID
						continue;
					}
				}
			} //末尾不换行
			if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
			{
				//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)
				if (boost::regex_search(mystr[8], cds_match_result, cds_match))
				{
					//printf("i am in cds...\n");
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3] << std::endl; //cds有name就用name
					out << mystr[0] << "\t" << mystr[1] << "\t" << "exon" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "." << "\t" << "ID=exon:" << cds_match_result[1] << "_rep" << ";Parent=" << cds_match_result[2] << std::endl; //cds有name就用name

					continue;
				}
				//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+
				if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
				{
					//printf("i am in cds1...\n");
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1] << std::endl; //没有Name的就等于cds的ID
					out << mystr[0] << "\t" << mystr[1] << "\t" << "exon" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "." << "\t" << "ID=exon:" << cds_match_result1[1] << "_rep" << ";Parent=" << cds_match_result1[2] << std::endl; //没有Name的就等于cds的ID
					continue;
				}
			}
			continue;
		}
		if (markmrna &&markexon) //没有gene和cds
		{
			if (mystr[2] == "mRNA" || mystr[2] == "mrna" || mystr[2] == "MRNA")
			{
				//([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)"; //这里没有parent，因为自己就是带头大哥
				//ID=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)
				if (boost::regex_search(mystr[8], result, pattern))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << "gene" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=gene:" << result[1] << "_rep" << ";Name=gene:" << result[1] << "_rep" << std::endl;
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result[1] << ";" << "Parent=gene:" << result[1] << "_rep" << ";Name=" << result[2] << ";gene=gene:" << result[1] << "_rep" << std::endl;
					continue;
				}
				//"ID=([^;\\s\"]++)[^\"]+" 没有Name
				if (boost::regex_search(mystr[8], gene_id, geneid))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << "gene" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=gene:" << gene_id[1] << "_rep" << ";Name=gene:" << gene_id[1] << "_rep" << std::endl; //name是自己的ID
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << gene_id[1] << ";Parent=gene:" << gene_id[1] << "_rep" << ";Name=" << gene_id[1] << ";gene =gene:" << gene_id[1] << "_rep" << std::endl; //name为自己的ID
					continue;
				}
				continue;
			}
			if (in.eof()) //末尾不换行
			{
				if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
				{
					//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+
					if (boost::regex_search(mystr[8], exon_result, exon_id))
					{
						if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 0)
						{
							out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "0" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
						}
						if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 1)
						{
							out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "1" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
						}
						if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 2)
						{
							out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "2" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
						}
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << exon_result[1] << ";Parent=" << exon_result[2];
						continue;
					}
				}
				if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
				{
					//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+
					if (boost::regex_search(mystr[8], exon_result, exon_id))
					{
						if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 0)
						{
							out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "0" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
						}
						if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 1)
						{
							out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "1" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
						}
						if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 2)
						{
							out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "2" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
						}
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << exon_result[1] << ";Parent=" << exon_result[2] << std::endl;
						continue;
					}
					continue;
				}

				continue;
			}
			if (markgene&&markexon) //没有mrna和cds
			{
				if (mystr[2] == "gene" || mystr[2] == "Gene" || mystr[2] == "GENE")
				{
					//"ID=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)", 
					if (boost::regex_search(mystr[8], result, pattern))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result[1] << ";" << "Name=" << result[2] << std::endl;
						out << mystr[0] << "\t" << mystr[1] << "\t" << "mRNA" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=mrna:" << result[1] << "_rep" << ";" << "Parent=" << result[1] << ";Name=" << "mrna:" << result[1] << "_rep" << ";gene=" << result[2] << std::endl;
						continue;
					}
					//ID=([^;\\s\"]++)[^\"]+  只有id没有name的gene
					if (boost::regex_search(mystr[8], gene_id, geneid))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << gene_id[1] << ";" << "Name=" << gene_id[1] << std::endl; //gene没有name的直接用id代替name
						out << mystr[0] << "\t" << mystr[1] << "\t" << "mRNA" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=mrna:" << result[1] << "_rep" << ";" << "Parent=" << result[1] << ";Name=" << "mrna:" << result[1] << "_rep" << ";gene=" << result[1] << std::endl; //gene没有name的直接用id代替name
					}
					continue;

				}
				if (in.eof()) //末尾不换行
				{
					if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
					{
						//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+
						if (boost::regex_search(mystr[8], exon_result, exon_id))
						{
							if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 0)
							{
								out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "0" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
							}
							if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 1)
							{
								out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "1" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
							}
							if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 2)
							{
								out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "2" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
							}
							out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << exon_result[1] << ";Parent=" << exon_result[2];
							continue;
						}
						continue;
					}
				}
				if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
				{
					//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+
					if (boost::regex_search(mystr[8], exon_result, exon_id))
					{
						if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 0)
						{
							out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "0" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
						}
						if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 1)
						{
							out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "1" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
						}
						if ((atoi(mystr[4].c_str()) - (atoi(mystr[3].c_str()) + 1)) % 3 == 2)
						{
							out << mystr[0] << "\t" << mystr[1] << "\t" << "CDS" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "2" << "\tID=cds:" << exon_result[1] << "_rep" << ";Parent=" << exon_result[2] << ";Name=cds:" << exon_result[1] << "_rep" << std::endl;
						}
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << exon_result[1] << ";Parent=" << exon_result[2] << std::endl;
						continue;
					}
					continue;
				}
				continue;
			}

		}
		in.close();
		out.close();
		return 1;
	}
}




int checkgff3update(const char *pathingff3, const char *pathinfa)
{
	using namespace boost;
	std::ifstream ingff3(pathingff3, std::ios::in);  //打开gff3文件
	std::ifstream infa(pathinfa, std::ios::in); //打开genome的fa文件
	std::map<std::string, std::string> genefa_seq; //用于存储fa文件
	std::string linefa; //用来存储fa流
	while (getline(infa, linefa, '>'))
	{
		std::string gene_id;
		std::string gene_seq;
		if (linefa.length() == 0 || linefa.find(0x0A) == 0)   //用来查找上一行是否是空或者回车,换行符都在0位置，0x0A=‘\n’
		{
			continue; //跳过开头的空行
		}
		size_t pos_n = linefa.find('\n'); //查找第一次换行符出现的位置，为了得到基因的ID;size_t是unsigned int类型
		gene_id.assign(linefa, 0, pos_n); //才是对的
		if (size_t pos_space = gene_id.find(0x20))  //找到space 空格
		{
			gene_id.assign(gene_id, 0, pos_space); //截取
		}
		if (size_t pos_tab = gene_id.find(0x09)) //找到水平制表符
		{
			gene_id.assign(gene_id, 0, pos_tab);
		}
		gene_seq.assign(linefa, pos_n + 1, linefa.length() - pos_n);  //获取序列
		// remove函数在头文件#include<algorithm>中，不包含这个头文件remove函数就不能使用
		//gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), '\n'), gene_seq.end()); 
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x0A), gene_seq.end()); //删除换行符，同上
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x09), gene_seq.end()); //删除“水平制表符”
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x08), gene_seq.end()); //删除backspace
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x20), gene_seq.end()); //删除空格
		transform(gene_seq.begin(), gene_seq.end(), gene_seq.begin(), ::toupper); //将字符串全部转换为大写 ，::tolower是小写
		std::cout << "gene_id " << gene_id << std::endl;
		//std::cout << "seq " << gene_seq << std::endl;
		genefa_seq[gene_id] = gene_seq;
	}
	infa.close(); //读取fa文件完毕
	//开始读取gff3
	std::string line0;
	std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,std::map<std::string,unsigned int>>> > > gene2trans; //保存基因和转录本对应信息 //尖括号必须分开写
	std::map<std::string, std::map<std::string, std::vector<std::string>> > trans_info; //保存exon和cds的位置信息   //尖括号必须分开写,vs编译器不用分开写，linux必须分开写
	//boost::regex pattern("ID=([^;\\s\"]{1,30});[^\\s\"]{0,50}Parent=([^;\\s\"]{1,30})", boost::regex::perl| boost::regex::icase); //少用.这样正则能匹配快很多,boost::regex::perl使用perl语法
	//++流模式只能放到捕获小括号里
	boost::regex pattern("ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)[^\\s\"]+gene=([\\w\\d]++);", boost::regex::perl | boost::regex::icase); //++是流模式减少回溯，提高正则速度
	boost::smatch result;
	boost::regex pattern1("Parent=([^;\\s\"]{1,30})", boost::regex::perl| boost::regex::icase);
	boost::smatch result1;
	boost::regex pattern2("Parent=([^;\\s\"]{1,30})", boost::regex::perl| boost::regex::icase); //大小写不敏感
	boost::smatch result2;
	//NC_003070.9	RefSeq	mRNA	3631	5899	.	+	.	ID=rna-NM_099983.2;Parent=gene-AT1G01010;Dbxref=Araport:AT1G01010,TAIR:AT1G01010,GeneID:839580,Genbank:NM_099983.2;Name=NM_099983.2;gbkey=mRNA;gene=NAC001;inference=similar to RNA sequence%2C mRNA:INSD:BT001115.1%2CINSD:AF439834.1%2CINSD:AK226863.1;locus_tag=AT1G01010;orig_protein_id=gnl|JCVI|AT1G01010.1;orig_transcript_id=gnl|JCVI|mRNA.AT1G01010.1;product=NAC domain containing protein 1;transcript_id=NM_099983.2
	while (getline(ingff3, line0))
	{
		if (line0.find("#") == 0) //跳过#开头的行
		{
			continue;
		}
		std::vector<std::string> mystr; //保存切割的字符串
		mystr = split(line0, "\t");
		std::cout << mystr[2] << std::endl;
		if (mystr[2] == "gene" || mystr[2] == "Gene" || mystr[2] == "GENE") //有些ggf3文件没有mRNA  大肠杆菌gff3文件就是这样
		{

			std::vector<std::string> mystr1 = split(mystr[8], ";"); //获取最后一行的分割字符串
			std::vector<std::string> myname = split(mystr1[2], "="); //获取name名
			std::vector<std::string> mygeneid = split(mystr1[0], "="); //获取基因名
		}
		if (mystr[2] == "mRNA" || mystr[2] == "mrna" || mystr[2] == "MRNA")
		{
			bool valid = boost::regex_search(mystr[8], result, pattern);
			//perl
			//$tmp[8] =~ /.*ID=([^;\s\"]+).*Parent=([^;\s\"]+)/;
			//std::string pattern(".*ID=([^;\\s\"]+).*Parent=([^;\\s\"]+);.*");
			if (valid)
			{
				std::cout << result[2] << std::endl;
				gene2trans[mystr[0]][result[2]][result[1]][result[3]][result[4]] = 1;
				trans_info[result[1]]["start"].push_back(mystr[3]);
				trans_info[result[1]]["end"].push_back(mystr[4]);
				trans_info[result[1]]["strand"].push_back(mystr[6]);
				std::cout << "mRNA " << result[1] << std::endl;
			}
		}
		if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
		{
			// $tmp[8] =~ /.*Parent=([^;\s\"]+)/;
			//std::string pattern1(".*Parent=([^;\\s\"]+)");
			bool valid = boost::regex_search(mystr[8], result1, pattern1);
			if (valid)
			{
				trans_info[result1[1]]["Exonarray"].push_back(mystr[3]);
				trans_info[result1[1]]["Exonarray"].push_back(mystr[4]); //result[0]是被匹配的字符串本身，不是捕获的字符串，1才是
				//std::cout << "result1[0]" << result1[1] << std::endl;
			}
		}
		if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
		{
			bool valid = boost::regex_search(mystr[8], result2, pattern2);
			//std::vector<std::string> temp2;
			// $tmp[8] =~ /.*Parent=([^;\s\"]+)/;
			//std::string pattern2(".*Parent=([^;\\s\"]+)");
			if (valid)
			{
				trans_info[result2[1]]["Cdsarray"].push_back(mystr[3]);
				trans_info[result2[1]]["Cdsarray"].push_back(mystr[4]);
				//std::cout<<mystr[3]<<std::endl;
				//std::cout << "result[0]" << result2[1] << std::endl;
			}
		}
	}
	ingff3.close(); //读取gff3文件完毕
	//输出文件
	std::string path1("C:\\Users\\xiaohui\\Desktop\\get_cds_seq\\gene2transcript.list");
	std::ofstream gene2transcript(path1,std::ios::out);
	std::cout << "triple keys" << std::endl;
	gene2transcript << "#Gene ID and Transcript ID of gff3 file statistic" << std::endl;
	gene2transcript << "#GeneID\t" << "TranscriptID" << "\tTranscript Name"<<"\tGene Name"<<std::endl;
	std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,std::map<std::string,unsigned int>>>>>::iterator chr;
	std::map<std::string, std::map<std::string, std::map<std::string,std::map<std::string,unsigned int>>>>::iterator gene;
	std::map<std::string, std::map<std::string,std::map<std::string,unsigned int>>>::iterator iso;
	std::map<std::string, std::map<std::string, unsigned int>>::iterator trans_name;
	std::map<std::string, unsigned int>::iterator gene_name;
	for (chr = gene2trans.begin(); chr != gene2trans.end(); chr++)
	{
		for (gene = chr->second.begin(); gene != chr->second.end(); gene++)
		{
			for (iso = gene->second.begin(); iso != gene->second.end(); iso++)
			{
				for (trans_name = iso->second.begin(); trans_name != iso->second.end(); trans_name++)
				{
					for (gene_name = trans_name->second.begin(); gene_name != trans_name->second.end(); gene_name++)
					{
						gene2transcript << gene->first << "\t" << iso->first <<"\t"<<trans_name->first<<"\t"<<gene_name->first<< std::endl;
					}
				}
				
			}
		}
	}
	gene2transcript.close();
	//
	std::string path2("C:\\Users\\xiaohui\\Desktop\\get_cds_seq\\Know.transcript.fa");
	std::ofstream transcript_fa(path2,std::ios::out);
	std::string path3("C:\\Users\\xiaohui\\Desktop\\get_cds_seq\\Know.longest_transcript.fa");
	std::ofstream longest_fa(path3, std::ios::out);
	std::string path4("C:\\Users\\xiaohui\\Desktop\\get_cds_seq\\Know.trans_max_id.list");
	std::ofstream trans_max_id1(path4,std::ios::out);
	std::string path5("C:\\Users\\xiaohui\\Desktop\\get_cds_seq\\Know.longest_cds.fa");
	std::ofstream cds_max_seq1(path5, std::ios::out);
	std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, unsigned int>>>>>::iterator chr1;
	std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, unsigned int>>>>::iterator gene1;
	std::map<std::string, std::map<std::string, std::map<std::string, unsigned int>>>::iterator iso1;
	std::map<std::string, std::map<std::string, unsigned int>>::iterator trans_name1;
	std::map<std::string, unsigned int>::iterator gene_name1;
	std::map<std::string, std::map<std::string, std::vector<std::string>> >::iterator ios1;
	//std::map<std::string, std::vector<std::string>>::iterator 
	for (chr1 = gene2trans.begin(); chr1 != gene2trans.end(); chr1++)
	{
		std::cout << "chromosome " << chr1->first << std::endl;
		for (gene1 = chr1->second.begin(); gene1 != chr1->second.end(); gene1++)
		{
			std::cout << "gene " << gene1->first << std::endl;
			std::string cds_max_seq;
			int trans_max_len = 0;
			std::string  trans_max_id;
			std::string trans_max_seq;
			int skip = 1;
			for (iso1 = gene1->second.begin(); iso1 != gene1->second.end(); iso1++)
			{
				std::cout << "iso1 " << iso1->first << std::endl;
				if (!(trans_info[iso1->first]["Exonarray"].size()>=1)) //如果存在这个键和值的话
				std::cout << "Number " << trans_info[iso1->first]["Cdsarray"][0] << std::endl;
				std::string trans_seq, cds_seq; int trans_len=0, cds_len=0;
				if (!(trans_info[iso1->first]["Exonarray"].size()>=1)) //如果存在这个键和值的话
				{
					std::cout<< "WARN: no exon in "<< iso1->first <<std::endl;
					continue;
				}
				std::vector<int> exon=sortstring(trans_info[iso1->first]["Exonarray"]); //排序从小到大
				std::vector<int> cds =sortstring(trans_info[iso1->first]["Cdsarray"]); //从小到大排序
				for (size_t i = 0; i < exon.size();i += 2) 
				{
					size_t len = (exon[i + 1] - exon[i] + 1);
					std::cout << "len " << len << std::endl;
					trans_len += len;
					std::cout << "chr1->first"<<chr1->first << std::endl;
					trans_seq = genefa_seq[chr1->first].substr(exon[i] - 1, len);
					std::cout << "seq " << genefa_seq[chr1->first] << std::endl;
				}
				for (size_t i = 0; i <cds.size(); i += 2)
				{
					size_t len = (cds[i + 1] - cds[i] + 1);
					cds_len += len;
					cds_seq += genefa_seq[chr1->first].substr(cds[i] - 1, len);
				}
				//std::vector<std::string>::iterator iter = std::find(trans_info[iso1->first]["strand"].begin(), trans_info[iso1->first]["strand"].end(), '-');//返回的是一个迭代器指针
				if (trans_info[iso1->first]["strand"][0]=="-") 
				{
					trans_seq=tr_reverse(trans_seq);
					//trans_seq = ~tr / ATCG / TAGC / ;
					//trans_seq = tr_reverse(trans_seq);
					//$trans_seq = ~tr / RYMKHDBV / YRKMDHVB / ;    //add by Liu Tao 2016 - 03 - 03, convert degenerate bases in raw sequence.
						//$trans_seq = reverse $trans_seq;
					cds_seq = tr_reverse(cds_seq);
					//$cds_seq = ~tr / RYMKHDBV / YRKMDHVB / ;dddddd
					//$cds_seq = reverse $cds_seq;
				}
				//if $trans_seq N% >= 80%
				int N_num = count(trans_seq.begin(), trans_seq.end(), 'N'); //统计序列中N的个数
				double N_ratio = (N_num*1.0) / (trans_seq.length())*1.0;
				if (N_ratio > 0.8)
				{
					std::cout<<"WARN: too much N in "<< iso1->first<<std::endl;
				}
				skip = 0;
				transcript_fa << ">" << ios1->first << "\n" << trans_seq << "\n";

				if (trans_max_len==0) {
					trans_max_len = trans_len;
					trans_max_id = iso1->first;
					trans_max_seq = trans_seq;

					cds_max_seq = cds_seq;
				}
				else if(trans_len > trans_max_len) {
					trans_max_len = trans_len;
					trans_max_id = iso1->first;
					trans_max_seq = trans_seq;
					cds_max_seq = cds_seq;
				}
				if (skip == 1) {
					continue;
				}
				//my $pep_max_seq = &TranslateDNASeq($cds_max_seq);
				longest_fa << ">" << gene1->first << "\n" << trans_max_seq << "\n";
				trans_max_id1 << gene1->first << "\t" << trans_max_id << "\n";
				cds_max_seq1 << gene1->first << "\n" << cds_max_seq << "\n";
				
				

			}
		}
	}
	longest_fa.close();
	trans_max_id1.close();
	cds_max_seq1.close();
	transcript_fa.close();
	return 1;
}




void gff3andgenome(const char *pathingff3, const char *pathingenome,const char *pathout)
{
	time_t start, end;
	time(&start); //获取开始的时间
	using std::endl;
	unsigned long count_del_slash_n = 0;
	unsigned long gene_count = 0; //统计基因个数
	unsigned long transcript_count = 0; //统计转录本数
	unsigned long exon_count = 0; //统计外显子数
	unsigned long cds_count = 0;//统计cds
	std::map<std::string, std::string> genomefa; //存储基因组fa文件
	std::ifstream ingff3(pathingff3, std::ios::in);  //打开gff3文件
	std::ifstream ingenome(pathingenome, std::ios::in); //打开genome文件
	//开始读取基因组文件，存入map容器
	std::string temp;
	while (getline(ingenome,temp,'>'))
	{
		if (temp.length() == 0 || temp.find(0x0A) == 0)   //用来查找上一行是否是空或者回车,换行符都在0位置，0x0A=‘\n’
		{
			continue; //跳过开头的空行
		}
		size_t pos_n = temp.find('\n'); //查找第一次换行符出现的位置，为了得到基因的ID;size_t是unsigned int类型
		std::string gene_id, gene_seq;
		gene_id.assign(temp, 0, pos_n); //才是对的,截取第一个换行前的字符串作为gene_id
		if (size_t pos_space = gene_id.find(0x20))  //找到space 空格
		{
			gene_id.assign(gene_id, 0, pos_space); //截取
		}
		if (size_t pos_tab = gene_id.find(0x09)) //找到水平制表符
		{
			gene_id.assign(gene_id, 0, pos_tab);
		}
		gene_seq.assign(temp, pos_n + 1, temp.length() - pos_n); //+1是为了去掉结尾的换行符号  //这个也可以用(获取序列)
		//remove函数在头文件#include<algorithm>中，不包含这个头文件remove函数就不能使用
		//gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), '\n'), gene_seq.end()); 
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x0A), gene_seq.end()); //删除换行符，同上
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x09), gene_seq.end()); //删除“水平制表符”
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x08), gene_seq.end()); //删除backspace
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x20), gene_seq.end()); //删除空格
		//genomefa[gene_id] = gene_seq;  //fa序列入容器map
		genomefa.insert(std::pair<std::string, std::string>(gene_id, gene_seq)); //这样也行
	}
	printf("fa读入内存完毕\n");
	//fa读入内存完毕

	std::ofstream out(pathout, std::ios::out);
	const char *pathlog = "C:\\Users\\xiaohui\\Desktop\\gff3.statistic";
	const char *path_longest_transcripts = "C:\\Users\\xiaohui\\Desktop\\Know.longest_transcripts.fa";
	std::ofstream longest_transcripts(path_longest_transcripts,std::ios::out);
	std::ofstream outlog(pathlog, std::ios::out);
	if (!ingff3.is_open() || !out.is_open())
	{
		printf("Your inputs files failed...\n");
		return;
	}
	outlog << "#gff3 statistic" << "\n";
	std::string line0;
	std::map<std::string, int> chrom;
	while (getline(ingff3, line0))
	{
		if (line0.find("#") == 0) //跳过#开头的行
		{
			continue;
		}
		std::vector<std::string> mystr;
		mystr = split(line0, "\t");
		if (mystr[2] == "gene" || mystr[2] == "Gene" || mystr[2] == "GENE")
		{
			chrom[mystr[0]]++; //统计染色体数目和一条染色体上有多少个基因
			gene_count++;
		}
		if (mystr[2] == "mRNA" || mystr[2] == "mrna" || mystr[2] == "MRNA")
		{
			transcript_count++;
		}
		if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
		{
			exon_count++;
		}
		if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
		{
			cds_count++;
		}
	}
	outlog << "Gene numbers:" << gene_count << "\n" << "Transcript numbers:" << transcript_count << "\n" << "Exon numbers:" << exon_count << "\n" << "CDS numbers:" << cds_count << endl;
	outlog << "Chromosome numbers:" << chrom.size() << endl; //打印染色体数目
	std::map<std::string, int>::iterator ibegin = chrom.begin();
	for (; ibegin != chrom.end(); ibegin++)
	{
		outlog << "Chromosome Name: " << ibegin->first << "\tGene numbers: " << ibegin->second << endl; //打印一条染色体上有多少条基因
	}
	ingff3.clear(std::ios::goodbit);//在seekg之前需要调用流对象的clear方法，把流的标记清除掉，清除以后就可以正常调用seekg方法了。
	ingff3.seekg(std::ios::beg); //begin的缩写,将文件指针重新回到头部 //seekp设置输出流文件指针位置，seekg设置输入流文件指针位置
	std::string line;
	std::string catcds;
	while (getline(ingff3, line)) //默认换行符为一行
	{
		if (line.find("#") == 0) //跳过#开头的行
		{
			continue;
		}
		std::vector<std::string> mystr;
		mystr = split(line, "\t");
		if (genomefa.count(mystr[0])>0)  //说明有这个键
		{
			if (mystr[2] == "gene" || mystr[2] == "Gene" || mystr[2] == "GENE")
			{
				std::vector<std::string> mystr1 = split(mystr[8], ";"); //获取最后一行的分割字符串
				std::vector<std::string> myname = split(mystr1[2], "="); //获取name名
				std::vector<std::string> mygeneid = split(mystr1[0], "="); //获取基因名
				catcds="";
				longest_transcripts << ">" << mygeneid[1] << " " << myname[1] << "\n";
			}
			if (mystr[2] == "mRNA" || mystr[2] == "mrna" || mystr[2] == "MRNA")
			{
				continue;
			}
			if (mystr[2] == "exon" || mystr[2] == "Exon" || mystr[2] == "EXON")
			{
				continue;
			}
			if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds") //将CDS拼接一起
			{
				unsigned int start = atoi(mystr[3].c_str());
				unsigned int end = atoi(mystr[4].c_str());
				catcds += genomefa[mystr[0]].substr(start,end-start+1);
			}
				longest_transcripts << catcds << endl;
			
		}
	}
		//"NC_000913.3	RefSeq	gene	190	255	.	+	.	ID=gene-b0001;Dbxref=ASAP:ABE-0000006,ECOCYC:EG11277,EcoGene:EG11277,GeneID:944742;Name=thrL;gbkey=Gene;gene=thrL;gene_biotype=protein_coding;gene_synonym=ECK0001;locus_tag=b0001"
	
	ingff3.close();
	out.close();
	outlog.close();
	longest_transcripts.close();
	time(&end); //获取结束的时间
	printf("Elapse %d min\n", (unsigned int)(end - start) / 60); //打印时间
}



