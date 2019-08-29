#define _CRT_SECURE_NO_WARNINGS
#include "function.h"


 inline std::vector<std::string> split(std::string str, std::string pattern)  //�ַ����и��������һ��vector���飨������������ǿ������ִ��Ч�ʣ�
{
		std::string::size_type pos;
		std::vector<std::string> result;
		str += pattern;//��չ�ַ����Է������
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
	return result; //����һ��vector����
}


 //get N number
 



void gff3statistic(const char *pathin,const char *pathout)
{
	time_t start, end;
	time(&start); //��ȡ��ʼ��ʱ��
	using std::endl;
	unsigned long count_del_slash_n = 0;
	unsigned long gene_count = 0; //ͳ�ƻ������
	unsigned long transcript_count = 0; //ͳ��ת¼����
	unsigned long exon_count = 0; //ͳ����������
	unsigned long cds_count = 0;//ͳ��cds
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
		if (line0.find("#") == 0) //����#��ͷ����
		{
			continue;
		}
		std::vector<std::string> mystr;
		mystr = split(line0, "\t");
		if (mystr[2] == "gene" || mystr[2] == "Gene" || mystr[2] == "GENE")
		{
			chrom[mystr[0]]++; //ͳ��Ⱦɫ����Ŀ��һ��Ⱦɫ�����ж��ٸ�����
			gene_count++;
		}
		if (mystr[2] == "mRNA" || mystr[2] == "mrna" || mystr[2] == "MRNA")
		{
			transnum[mystr[0]]++; //ͳ��һ��Ⱦɫ���϶��ٸ�ת¼��
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
	outlog << "Chromosome numbers:" << chrom.size() << endl; //��ӡȾɫ����Ŀ
	std::map<std::string, int>::iterator ibegin = chrom.begin();
	std::map<std::string, unsigned int>::iterator itranbegin = transnum.begin();
	for (; ibegin != chrom.end(); ++ibegin)
	{
			outlog << "Chromosome Name: " << ibegin->first << "\tGene numbers: " << ibegin->second << "\t"<<transnum[ibegin->first]<<endl;// "\t" << itranbegin->second << endl; //��ӡһ��Ⱦɫ�����ж���������Ͷ��ٸ�ת¼��
		//outlog << "Chromosome Name: " << ibegin->first << "\tGene numbers: " << ibegin->second << endl;// "\t" << itranbegin->second << endl; //��ӡһ��Ⱦɫ�����ж���������Ͷ��ٸ�ת¼��
			
	}
	in.clear(std::ios::goodbit);//��seekg֮ǰ��Ҫ�����������clear�����������ı�������������Ժ�Ϳ�����������seekg�����ˡ�
	in.seekg(std::ios::beg); //begin����д,���ļ�ָ�����»ص�ͷ�� //seekp����������ļ�ָ��λ�ã�seekg�����������ļ�ָ��λ��
	std::string line;
	while (getline(in, line)) //Ĭ�ϻ��з�Ϊһ��
	{
		if (line.find("#") == 0) //����#��ͷ����
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
			out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << mystr1[0] << ":" << myname[1] << ";" << mystr1[1]<< endl; //mystr1[1]��parent
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
				out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << mystr1[0]<< ";" << mystr1[1]<< endl; //������û��Name
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
	time(&end); //��ȡ������ʱ��
	printf("Elapse %d min\n", (unsigned int)(end - start) / 60); //��ӡʱ��
}

//������ʽ����
inline bool is_regex(const std::string& object,std::vector<std::string> &myresult, const std::string &str)
{
	using std::cout;
	using std::endl;
	const std::regex pattern(str);
	std::smatch result;
	bool valid = std::regex_search(object, result, pattern);
	//�˴�result�������п��ޣ�result��һ���ַ������飬�����洢������ʽ�������ŵ����ݡ�
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
		//for (size_t i = 1; i < result.size(); i++) //����Ŀ���ַ���
		//{
		//	cout << result[i]<<endl ;
		//	myresult.push_back(result[i]); //��һ����Ŀ���ַ��������ڶ����ǵ�һ��С���ţ���������
	//	}
	
	return valid;
}

//string �ַ���ת��Ϊ���ֺ�����
std::vector<int> sortstring(std::vector<std::string> &in)
{
	std::vector<int> myreturn;
	printf("i am in\n");
	for (size_t i = 0; i < in.size(); i++)
	{
		myreturn.push_back(atoi(in[i].c_str()));
	}
	std::sort(myreturn.begin(), myreturn.end(), std::less<int>()); //��������
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
int check_no_coding_gene(const char *ingff3,size_t &linenum)   //linenum����ͳ���ļ�����
{
	enum {correct=1,incorrect=0,rRNA1=2,tRNA1=3,lncRNA1=4,miRN1A=5,circRNA1=6,noexon=7};  //��������ö������
	short markgene=0, markmrna=0, markexon=0, markcds=0,rRNA=0,tRNA=0,LncRNA=0,CircRNA=0,MiRNA=0;
	int rRNA_count = 0, tRNA_count = 0, circRNA_count = 0, LncRNA_count = 0, MiRNA_count = 0, cds_count = 0, gene_count = 0,mRNA_count=0;
	short no_coding_rna_exon_count= 0;
	std::map<std::string, size_t>stat_exon;
	std::map<std::string, size_t> stat_rna;
	std::map<std::string, size_t> stat_gene;
	//rRNAƥ��
	//boost::regex rRNA_pattern("R[ibosomal\\s]{0,9}RNA", boost::regex::perl | boost::regex::icase);
	//boost::smatch rRNA_result; //���ñ�����Ľ�����Ͳ��ö����������
	boost::smatch markresult;
	std::string xrna;
	std::ifstream in(ingff3, std::ios::in);
	if (!in.is_open())
	{
		std::cout << "Open file *.gff3 failed." << std::endl;
		return incorrect; //����һ��ö�����ͣ�����ĳ���ֵ��0��defaultΪint����
	}
	std::string line;
	while (getline(in, line))
	{
		if (line.find("#") == 0 || line.length() == 0 || line.find(0x0A) == 0) //����#��ͷ���� //����������һ���Ƿ��ǿջ��߻س�,���з�����0λ�ã�0x0A=��\n��
		{
			continue;
		}
		std::vector<std::string> mystr(split(line,"\t")); //�����и���ַ���
		if (boost::regex_match(mystr[2], markresult, boost::regex("(gene)", boost::regex::perl | boost::regex::icase)))
		{
			gene_count++;
			markgene = 1;
			continue;

		}
		//xrna_result Ҫ���ȶ���
		if (boost::regex_match(mystr[2], markresult,boost::regex("([\\w_-\\s]{0,15}RNA)", boost::regex::perl | boost::regex::icase))) //ƥ�䵽*RNA,����rRNA,tRNA��LncRNA�ȵ�
		{
			markmrna = 1;
			stat_rna[markresult[1]]++; //ͳ�Ʋ�ͬ��mRNA����
		}
		if (boost::regex_match(mystr[2], markresult, boost::regex("(exon)", boost::regex::perl | boost::regex::icase))) //ƥ��exon
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
		//match����ȫƥ�䣬search�������Ӵ�
		if (boost::regex_match(mystr[2],boost::regex("R[ibosomal_-\\s]{0,10}RNA", boost::regex::perl | boost::regex::icase))) //������RNA��ֻ��exonû��cds ,��Ҫ��ƥ��Ľ���ŵ�һ��������Ͳ��õڶ�������result
		{
			xrna = mystr[2];
			rRNA_count++;
			rRNA = 1;
			continue;
		}
		if (boost::regex_match(mystr[2], boost::regex("T[ransfer\\s_-ribonucleic_-\\sacid_-\\sRNA]{0,34}", boost::regex::perl | boost::regex::icase))) //ת��RNAֻ��exonû��cds
		{
			xrna = mystr[2];
			tRNA_count++;
			tRNA = 1;
			continue;
		}
		if (boost::regex_match(mystr[2], boost::regex("l[ong_-\\snoncoding_-\\s]{0,19}RNA", boost::regex::perl | boost::regex::icase))) //�����Ǳ���
		{
			xrna = mystr[2];
			LncRNA_count++;
			LncRNA = 1;
			continue;
		}
		if (boost::regex_match(mystr[2], boost::regex("Mi[cro_-\\s]{0,9}RNA", boost::regex::perl | boost::regex::icase))) //СRNA
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
			stat_exon[xrna]++;  //ͳ�Ʋ�ͬrna��exon����
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
	if (!(markexon==1)) //���û�������ӣ�Ҳ����ֻ��cds
	{
		std::cout << "gff3 have no exon mark.\n" << std::endl;
		if (markgene &&markmrna &&markcds) //��gene��mrna��cds
		{
			linenum = cds_count + mRNA_count + gene_count;
			return noexon;
		}
		if (markgene &&markcds) //��gene��cds
		{
			linenum = cds_count +gene_count;
			return noexon;
		}
		if (markmrna &&markcds) //��mrna��cds
		{
			linenum = cds_count + mRNA_count;
			return noexon;
		}
	}

	//�����������
	if (markgene && markmrna && markexon && markcds &&CircRNA &&MiRNA &&LncRNA &&rRNA &&tRNA) //��Щȫ������
	{
		std::map<std::string, size_t>::iterator begin = stat_exon.begin();
		for (; begin != stat_exon.end(); begin++)
		{
			linenum += begin->second;
		}
		linenum += cds_count + circRNA_count + MiRNA_count + LncRNA_count + rRNA_count + tRNA_count + mRNA_count + gene_count; //���û�еľ���0
		return 1;
	}
	//
	if (markgene && markmrna && markexon && markcds &&CircRNA &&MiRNA &&rRNA&&tRNA) //ֻ��crna,rrna,mirna,rRNA
	{
		std::map<std::string, size_t>::iterator begin = stat_exon.begin();
		for (; begin != stat_exon.end(); begin++)
		{
			linenum += begin->second;
		}
		linenum += cds_count + circRNA_count + MiRNA_count +rRNA_count + tRNA_count + mRNA_count + gene_count;
		return 1;
	}
	if (markgene && markmrna && markexon && markcds &&CircRNA &&MiRNA &&rRNA&&tRNA) //ֻ��crna,rrna,mirna,rRNA
	{
		std::map<std::string, size_t>::iterator begin = stat_exon.begin();
		for (; begin != stat_exon.end(); begin++)
		{
			linenum += begin->second;
		}
		linenum += cds_count + circRNA_count + MiRNA_count + rRNA_count + tRNA_count + mRNA_count + gene_count;
		return 1;
	}
	if (markgene && markmrna && markexon && markcds &&CircRNA) //ֻ��circrna
	{
		std::map<std::string, size_t>::iterator begin = stat_exon.begin();
		for (; begin != stat_exon.end(); begin++)
		{
			linenum += begin->second;
		}
		linenum += cds_count + circRNA_count + mRNA_count + gene_count; //���û�еľ���0
		return 1;
	}
	if (markgene && markmrna && markexon && markcds &&MiRNA) //ֻ��MiRNA
	{
		std::map<std::string, size_t>::iterator begin = stat_exon.begin();
		for (; begin != stat_exon.end(); begin++)
		{
			linenum += begin->second;
		}
		linenum += cds_count + MiRNA_count + mRNA_count + gene_count; //���û�еľ���0
		return 1;
	}
	if (markgene && markmrna && markexon && markcds &&rRNA) //ֻ��rrna
	{
		std::map<std::string, size_t>::iterator begin = stat_exon.begin();
		for (; begin != stat_exon.end(); begin++)
		{
			linenum += begin->second;
		}
		linenum += cds_count + rRNA_count + mRNA_count + gene_count; //���û�еľ���0
		return 1;
	}
	if (markgene && markmrna && markexon && markcds && LncRNA) //ֻ��lncrna
	{
		std::map<std::string, size_t>::iterator begin = stat_exon.begin();
		for (; begin != stat_exon.end(); begin++)
		{
			linenum += begin->second;
		}
		linenum += cds_count + LncRNA_count + mRNA_count + gene_count; //���û�еľ���0
		return 1;
	}
	if (markgene && markmrna && markexon && markcds && tRNA) //ֻ��trna
	{
		std::map<std::string, size_t>::iterator begin = stat_exon.begin();
		for (; begin != stat_exon.end(); begin++)
		{
			linenum += begin->second;
		}
		linenum += cds_count + tRNA_count + mRNA_count + gene_count; //���û�еľ���0
		return 1;
	}
	//û��gene
	if (markmrna && markexon && markcds &&CircRNA &&MiRNA &&LncRNA &&rRNA &&tRNA)
	{

	}
	if (markmrna && markexon && markcds && tRNA) //ת��rna
	{

	}
}









int checkgff3mark(const char *pathingff3, const char *outgff3)  //����һ����������������ʽ��gff3�ļ���Ϊ���溯��������
{
	int markgene = 0, flaggene = 0;
	int markmrna = 0, flagmrna = 0;
	int markexon = 0;
	int markcds = 0;
	int count = 1000; //��ȡ1000���ж��ļ��Ƿ��Ǳ�׼��gff3�ļ�
	std::ifstream in(pathingff3, std::ios::in);
	std::string line;
	boost::regex pattern("ID=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)", boost::regex::perl | boost::regex::icase); //++����ģʽ���ٻ��ݣ���������ٶ�
	boost::smatch result;
	boost::regex pattern1("ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)[^\\s\"]+gene=([\\w\\d]++);", boost::regex::perl | boost::regex::icase); //++����ģʽ���ٻ��ݣ���������ٶ�
	boost::smatch result1;
	boost::regex pattern2("ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)", boost::regex::perl | boost::regex::icase);
	boost::smatch result2;
	//mrna ��ƥ��
	boost::regex pattern_mrna("ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)", boost::regex::perl | boost::regex::icase);
	boost::smatch result_mrna;
	boost::regex pattern_mrna1("ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+", boost::regex::perl | boost::regex::icase); //û��name��ƥ��
	boost::smatch result_mrna1;
	//gene only id
	boost::regex geneid("ID=([^;\\s\"]++)[^\"]+", boost::regex::perl | boost::regex::icase);
	boost::smatch gene_id;
	//exon ƥ��
	boost::regex exon_id("ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]++)", boost::regex::perl | boost::regex::icase);
	boost::smatch exon_result;
	//ƥ��cds
	boost::regex cds_match("ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)", boost::regex::perl | boost::regex::icase);
	boost::smatch cds_match_result;
	boost::regex cds_match1("ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+", boost::regex::perl | boost::regex::icase);
	boost::smatch cds_match_result1;
	//û��ID��������,Name��������ID
	boost::regex exon_no_id("Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)", boost::regex::perl | boost::regex::icase);
	boost::smatch exon_no_id_result;
	std::map<std::string, std::string> gene;
	std::map<std::string, std::map<std::string, std::string>> mrnaname;
	std::map< std::string, std::map<std::string, std::string >> cdsname;
	std::cout << "markgene=" << markgene << "\tmarkmrna=" << markmrna << "\tmarkexon=" << markexon << "\tmarkcds=" << markcds << std::endl;

	while (getline(in, line))
	{

		if (line.find("#") == 0 || line.length() == 0 || line.find(0x0A) == 0) //����#��ͷ���� //����������һ���Ƿ��ǿջ��߻س�,���з�����0λ�ã�0x0A=��\n��
		{
			continue;
		}
		if (count == 0)
		{
			break; //����ѭ��
		}
		count--; //ѭ��50��//1000
		std::vector<std::string> mystr(split(line, "\t")); //�����и���ַ���
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
		if (mystr[2] == "rRNA" || mystr[2] == "rrna" || mystr[2] == "RRNA") //������RNA��ֻ��exonû��cds
		{

		}
		if (mystr[2] == "tRNA" || mystr[2] == "trna" || mystr[2] == "TRNA") //ת��RNAֻ��exonû��cds
		{

		}
		if (mystr[2] == "LncRNA" || mystr[2] == "lncrna" || mystr[2] == "LNCRNA") //�����Ǳ���
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
	in.clear(std::ios::goodbit);//��seekg֮ǰ��Ҫ�����������clear�����������ı�������������Ժ�Ϳ�����������seekg�����ˡ�
	in.seekg(std::ios::beg); //begin����д,���ļ�ָ�����»ص�ͷ�� //seekp����������ļ�ָ��λ�ã�seekg�����������ļ�ָ��λ��
	std::ofstream out(outgff3, std::ios::out);
	std::string line1;
	int flag1 = 0;
	std::string SaveGenenName;
	std::string SaveGeneID;
	while (getline(in, line1))
	{
		if (line1.find("#") == 0 || line1.length() == 0 || line1.find(0x0A) == 0) //����#��ͷ���� //����������һ���Ƿ��ǿջ��߻س�,���з�����0λ�ã�0x0A=��\n��
		{
			continue;
		}
		std::vector<std::string> mystr(split(line1, "\t")); //�����и���ַ���
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
				//ID=([^;\\s\"]++)[^\"]+  ֻ��idû��name��gene
				if (boost::regex_search(mystr[8], gene_id, geneid))
				{
					//printf("i am in third...\n");
					SaveGeneID = gene_id[1];
					flag1 = 2;
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << gene_id[1] << ";" << "Name=" << gene_id[1] << std::endl; //geneû��name��ֱ����id����name
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
				//"ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)" û��Name
				if (boost::regex_search(mystr[8], result2, pattern2))
				{
					//std::cout <<"result2="<< result2[0] << std::endl;
					//std::cout << "result2="<<result2[1] << std::endl;
					//std::cout << "result2="<<result2[2] << std::endl;
					//std::cout << "result="<<result[2] << std::endl;
					if (flag1 == 1)
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result2[1] << ";Parent=" << result2[2] << ";Name=" << result2[1] << ";gene=" << SaveGenenName << std::endl; //��name�͵���gene name
						flag1 = 0;
						continue;
					}
					if (flag1 == 2)
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result2[1] << ";Parent=" << result2[2] << ";Name=" << result2[1] << ";gene=" << SaveGeneID << std::endl; //ûname�͵���gene id
						flag1 = 0;
						continue;
					}
				}
				continue;
			}
			if (in.eof()) //�ж�����Ƕ�ȡ�����һ�У���ô���벻��ӡ����
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
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3]; //cds��name����name
						continue;
					}
					if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1]; //û��Name�ľ͵���cds��ID
						continue;
					}
				}
			} //��������һ�в���ӡ���з�

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
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3] << std::endl; //cds��name����name
					continue;
				}
				if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1] << std::endl; //û��Name�ľ͵���cds��ID
					continue;
				}
			}
			continue; //���ִ���������������������
		}
		if (markgene && markcds &&markexon) //��gene,cds��exon��û��mrnaִ��������
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
				//ID=([^;\\s\"]++)[^\"]+  ֻ��idû��name��gene
				if (boost::regex_search(mystr[8], gene_id, geneid))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << gene_id[1] << ";" << "Name=" << gene_id[1] << std::endl; //geneû��name��ֱ����id����name
					out << mystr[0] << "\t" << mystr[1] << "\t" << "mRNA" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=mrna:" << result[1] << "_rep" << ";" << "Parent=" << result[1] << ";Name=" << "mrna:" << result[1] << "_rep" << ";gene=" << result[1] << std::endl; //geneû��name��ֱ����id����name
				}
				continue;
			}
			if (in.eof()) //�ж�����Ƕ�ȡ�����һ�У���ô���벻��ӡ����
			{
				if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
				{
					if (boost::regex_search(mystr[8], cds_match_result, cds_match))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3]; //cds��name����name
						continue;
					}
					if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1]; //û��Name�ľ͵���cds��ID
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
			} //�����ĩβ�ˣ�����ӻ���

			if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
			{
				if (boost::regex_search(mystr[8], cds_match_result, cds_match))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3] << std::endl; //cds��name����name
					continue;
				}
				if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1] << std::endl; //û��Name�ľ͵���cds��ID
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
		if (markgene &&markmrna&& markcds) //��gene,mrna��exon��û��cdsִ�����
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
				//ID=([^;\\s\"]++)[^\"]+  ֻ��idû��name��gene
				if (boost::regex_search(mystr[8], gene_id, geneid))
				{
					flag1 = 2;
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << gene_id[1] << ";" << "Name=" << gene_id[1] << std::endl; //geneû��name��ֱ����id����name
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
				//"ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)" û��Name
				if (boost::regex_search(mystr[8], result2, pattern2))
				{
					if (flag1 == 1)
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result2[1] << ";Parent=" << result2[2] << ";Name=" << result2[1] << ";gene =" << result[2] << std::endl; //��name�͵���gene name
						flag1 = 0;
						continue;
					}
					if (flag1 == 2)
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result2[1] << ";Parent=" << result2[2] << ";Name=" << result2[1] << ";gene =" << result[1] << std::endl; //ûname�͵���gene id
						flag1 = 0;
						continue;
					}
				}
				continue;
			}
			if (in.eof()) //��������ļ�ĩβ�ˣ�������
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
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << exon_result[1] << ";Parent=" << exon_result[2]; //������
						continue;
					}
					continue;
				}
			}//������

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
		if (markmrna &&markcds &&markexon) //��mrna,exon��cdsû��gene
		{
			if (mystr[2] == "mRNA" || mystr[2] == "mrna" || mystr[2] == "MRNA")
			{
				//([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)"; //����û��parent����Ϊ�Լ����Ǵ�ͷ���
				//ID=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)
				if (boost::regex_search(mystr[8], result, pattern))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << "gene" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=gene:" << result[1] << "_rep" << ";Name=gene:" << result[1] << "_rep" << std::endl;
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result[1] << ";" << "Parent=gene:" << result[1] << "_rep" << ";Name=" << result[2] << ";gene=gene:" << result[1] << "_rep" << std::endl;
					continue;
				}
				//"ID=([^;\\s\"]++)[^\"]+" û��Name
				if (boost::regex_search(mystr[8], gene_id, geneid))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << "gene" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=gene:" << gene_id[1] << "_rep" << ";Name=gene:" << gene_id[1] << "_rep" << std::endl; //name���Լ���ID
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << gene_id[1] << ";Parent=gene:" << gene_id[1] << "_rep" << ";Name=" << gene_id[1] << ";gene =gene:" << gene_id[1] << "_rep" << std::endl; //nameΪ�Լ���ID
					continue;
				}
				continue;
			}
			if (in.eof()) //��������ļ�ĩβ�ˣ�������
			{
				if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
				{
					if (boost::regex_search(mystr[8], cds_match_result, cds_match))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3]; //cds��name����name
						continue;
					}
					if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1]; //û��Name�ľ͵���cds��ID
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
			} //������
			if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
			{
				if (boost::regex_search(mystr[8], cds_match_result, cds_match))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3] << std::endl; //cds��name����name
					continue;
				}
				if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1] << std::endl; //û��Name�ľ͵���cds��ID
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
		//������
		if (markmrna &&markcds) //��mrna��cdsû��gene��exon
		{
			if (mystr[2] == "mRNA" || mystr[2] == "mrna" || mystr[2] == "MRNA")
			{
				//([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)"; //����û��parent����Ϊ�Լ����Ǵ�ͷ���
				//ID=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)
				if (boost::regex_search(mystr[8], result, pattern))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << "gene" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=gene:" << result[1] << "_rep" << ";Name=gene:" << result[1] << "_rep" << std::endl;
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result[1] << ";" << "Parent=gene:" << result[1] << "_rep" << ";Name=" << result[2] << ";gene=gene:" << result[1] << "_rep" << std::endl;
					continue;
				}
				//"ID=([^;\\s\"]++)[^\"]+" û��Name
				if (boost::regex_search(mystr[8], gene_id, geneid))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << "gene" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=gene:" << gene_id[1] << "_rep" << ";Name=gene:" << gene_id[1] << "_rep" << std::endl; //name���Լ���ID
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << gene_id[1] << ";Parent=gene:" << gene_id[1] << "_rep" << ";Name=" << gene_id[1] << ";gene =gene:" << gene_id[1] << "_rep" << std::endl; //nameΪ�Լ���ID
					continue;
				}
				continue;
			}
			if (in.eof()) //��������ļ�ĩβ�ˣ�������
			{
				if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
				{
					//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)
					if (boost::regex_search(mystr[8], cds_match_result, cds_match))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3] << std::endl; //cds��name����name
						out << mystr[0] << "\t" << mystr[1] << "\t" << "exon" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "." << "\tID=exon:" << cds_match_result[1] << "_rep" << ";Parent=" << cds_match_result[2]; //cds��name����name

						continue;
					}
					//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+
					if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1] << std::endl; //û��Name�ľ͵���cds��ID
						out << mystr[0] << "\t" << mystr[1] << "\t" << "exon" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "." << "\tID=exon:" << cds_match_result1[1] << "_rep" << ";Parent=" << cds_match_result1[2]; //û��Name�ľ͵���cds��ID
						continue;
					}
				}
			} //ĩβ������

			if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
			{
				//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)
				if (boost::regex_search(mystr[8], cds_match_result, cds_match))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3] << std::endl; //cds��name����name
					out << mystr[0] << "\t" << mystr[1] << "\t" << "exon" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "." << "\tID=exon:" << cds_match_result[1] << "_rep" << ";Parent=" << cds_match_result[2] << std::endl; //cds��name����name

					continue;
				}
				//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+
				if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\tID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1] << std::endl; //û��Name�ľ͵���cds��ID
					out << mystr[0] << "\t" << mystr[1] << "\t" << "exon" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "." << "\tID=exon:" << cds_match_result1[1] << "_rep" << ";Parent=" << cds_match_result1[2] << std::endl; //û��Name�ľ͵���cds��ID
					continue;
				}
			}

			continue;
		}
		if (markgene &&markcds) //��gene��cdsû��mrna��exon
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
				//ID=([^;\\s\"]++)[^\"]+  ֻ��idû��name��gene
				if (boost::regex_search(mystr[8], gene_id, geneid))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << gene_id[1] << ";" << "Name=" << gene_id[1] << std::endl; //geneû��name��ֱ����id����name
					out << mystr[0] << "\t" << mystr[1] << "\t" << "mRNA" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=mrna:" << result[1] << "_rep" << ";" << "Parent=" << result[1] << ";Name=" << "mrna:" << result[1] << "_rep" << ";gene=" << result[1] << std::endl; //geneû��name��ֱ����id����name
				}
				continue;

			}
			if (in.eof()) //ĩβ������
			{
				if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
				{
					//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)
					if (boost::regex_search(mystr[8], cds_match_result, cds_match))
					{
						//printf("i am in cds...\n");
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3] << std::endl; //cds��name����name
						out << mystr[0] << "\t" << mystr[1] << "\t" << "exon" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "." << "\t" << "ID=exon:" << cds_match_result[1] << "_rep" << ";Parent=" << cds_match_result[2]; //cds��name����name

						continue;
					}
					//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+
					if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
					{
						//printf("i am in cds1...\n");
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1] << std::endl; //û��Name�ľ͵���cds��ID
						out << mystr[0] << "\t" << mystr[1] << "\t" << "exon" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "." << "\t" << "ID=exon:" << cds_match_result1[1] << "_rep" << ";Parent=" << cds_match_result1[2]; //û��Name�ľ͵���cds��ID
						continue;
					}
				}
			} //ĩβ������
			if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds")
			{
				//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)
				if (boost::regex_search(mystr[8], cds_match_result, cds_match))
				{
					//printf("i am in cds...\n");
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << cds_match_result[1] << ";Parent=" << cds_match_result[2] << ";Name=" << cds_match_result[3] << std::endl; //cds��name����name
					out << mystr[0] << "\t" << mystr[1] << "\t" << "exon" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "." << "\t" << "ID=exon:" << cds_match_result[1] << "_rep" << ";Parent=" << cds_match_result[2] << std::endl; //cds��name����name

					continue;
				}
				//ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+
				if (boost::regex_search(mystr[8], cds_match_result1, cds_match1))
				{
					//printf("i am in cds1...\n");
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << cds_match_result1[1] << ";Parent=" << cds_match_result1[2] << ";Name=" << cds_match_result1[1] << std::endl; //û��Name�ľ͵���cds��ID
					out << mystr[0] << "\t" << mystr[1] << "\t" << "exon" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << "." << "\t" << "ID=exon:" << cds_match_result1[1] << "_rep" << ";Parent=" << cds_match_result1[2] << std::endl; //û��Name�ľ͵���cds��ID
					continue;
				}
			}
			continue;
		}
		if (markmrna &&markexon) //û��gene��cds
		{
			if (mystr[2] == "mRNA" || mystr[2] == "mrna" || mystr[2] == "MRNA")
			{
				//([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]+)[^\\s\"]+Name=([^;\\s\"]++)"; //����û��parent����Ϊ�Լ����Ǵ�ͷ���
				//ID=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)
				if (boost::regex_search(mystr[8], result, pattern))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << "gene" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=gene:" << result[1] << "_rep" << ";Name=gene:" << result[1] << "_rep" << std::endl;
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << result[1] << ";" << "Parent=gene:" << result[1] << "_rep" << ";Name=" << result[2] << ";gene=gene:" << result[1] << "_rep" << std::endl;
					continue;
				}
				//"ID=([^;\\s\"]++)[^\"]+" û��Name
				if (boost::regex_search(mystr[8], gene_id, geneid))
				{
					out << mystr[0] << "\t" << mystr[1] << "\t" << "gene" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=gene:" << gene_id[1] << "_rep" << ";Name=gene:" << gene_id[1] << "_rep" << std::endl; //name���Լ���ID
					out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << gene_id[1] << ";Parent=gene:" << gene_id[1] << "_rep" << ";Name=" << gene_id[1] << ";gene =gene:" << gene_id[1] << "_rep" << std::endl; //nameΪ�Լ���ID
					continue;
				}
				continue;
			}
			if (in.eof()) //ĩβ������
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
			if (markgene&&markexon) //û��mrna��cds
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
					//ID=([^;\\s\"]++)[^\"]+  ֻ��idû��name��gene
					if (boost::regex_search(mystr[8], gene_id, geneid))
					{
						out << mystr[0] << "\t" << mystr[1] << "\t" << mystr[2] << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=" << gene_id[1] << ";" << "Name=" << gene_id[1] << std::endl; //geneû��name��ֱ����id����name
						out << mystr[0] << "\t" << mystr[1] << "\t" << "mRNA" << "\t" << mystr[3] << "\t" << mystr[4] << "\t" << mystr[5] << "\t" << mystr[6] << "\t" << mystr[7] << "\t" << "ID=mrna:" << result[1] << "_rep" << ";" << "Parent=" << result[1] << ";Name=" << "mrna:" << result[1] << "_rep" << ";gene=" << result[1] << std::endl; //geneû��name��ֱ����id����name
					}
					continue;

				}
				if (in.eof()) //ĩβ������
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
	std::ifstream ingff3(pathingff3, std::ios::in);  //��gff3�ļ�
	std::ifstream infa(pathinfa, std::ios::in); //��genome��fa�ļ�
	std::map<std::string, std::string> genefa_seq; //���ڴ洢fa�ļ�
	std::string linefa; //�����洢fa��
	while (getline(infa, linefa, '>'))
	{
		std::string gene_id;
		std::string gene_seq;
		if (linefa.length() == 0 || linefa.find(0x0A) == 0)   //����������һ���Ƿ��ǿջ��߻س�,���з�����0λ�ã�0x0A=��\n��
		{
			continue; //������ͷ�Ŀ���
		}
		size_t pos_n = linefa.find('\n'); //���ҵ�һ�λ��з����ֵ�λ�ã�Ϊ�˵õ������ID;size_t��unsigned int����
		gene_id.assign(linefa, 0, pos_n); //���ǶԵ�
		if (size_t pos_space = gene_id.find(0x20))  //�ҵ�space �ո�
		{
			gene_id.assign(gene_id, 0, pos_space); //��ȡ
		}
		if (size_t pos_tab = gene_id.find(0x09)) //�ҵ�ˮƽ�Ʊ��
		{
			gene_id.assign(gene_id, 0, pos_tab);
		}
		gene_seq.assign(linefa, pos_n + 1, linefa.length() - pos_n);  //��ȡ����
		// remove������ͷ�ļ�#include<algorithm>�У����������ͷ�ļ�remove�����Ͳ���ʹ��
		//gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), '\n'), gene_seq.end()); 
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x0A), gene_seq.end()); //ɾ�����з���ͬ��
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x09), gene_seq.end()); //ɾ����ˮƽ�Ʊ����
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x08), gene_seq.end()); //ɾ��backspace
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x20), gene_seq.end()); //ɾ���ո�
		transform(gene_seq.begin(), gene_seq.end(), gene_seq.begin(), ::toupper); //���ַ���ȫ��ת��Ϊ��д ��::tolower��Сд
		std::cout << "gene_id " << gene_id << std::endl;
		//std::cout << "seq " << gene_seq << std::endl;
		genefa_seq[gene_id] = gene_seq;
	}
	infa.close(); //��ȡfa�ļ����
	//��ʼ��ȡgff3
	std::string line0;
	std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,std::map<std::string,unsigned int>>> > > gene2trans; //��������ת¼����Ӧ��Ϣ //�����ű���ֿ�д
	std::map<std::string, std::map<std::string, std::vector<std::string>> > trans_info; //����exon��cds��λ����Ϣ   //�����ű���ֿ�д,vs���������÷ֿ�д��linux����ֿ�д
	//boost::regex pattern("ID=([^;\\s\"]{1,30});[^\\s\"]{0,50}Parent=([^;\\s\"]{1,30})", boost::regex::perl| boost::regex::icase); //����.����������ƥ���ܶ�,boost::regex::perlʹ��perl�﷨
	//++��ģʽֻ�ܷŵ�����С������
	boost::regex pattern("ID=([^;\\s\"]++)[^\\s\"]+Parent=([^;\\s\"]++)[^\\s\"]+Name=([^;\\s\"]++)[^\\s\"]+gene=([\\w\\d]++);", boost::regex::perl | boost::regex::icase); //++����ģʽ���ٻ��ݣ���������ٶ�
	boost::smatch result;
	boost::regex pattern1("Parent=([^;\\s\"]{1,30})", boost::regex::perl| boost::regex::icase);
	boost::smatch result1;
	boost::regex pattern2("Parent=([^;\\s\"]{1,30})", boost::regex::perl| boost::regex::icase); //��Сд������
	boost::smatch result2;
	//NC_003070.9	RefSeq	mRNA	3631	5899	.	+	.	ID=rna-NM_099983.2;Parent=gene-AT1G01010;Dbxref=Araport:AT1G01010,TAIR:AT1G01010,GeneID:839580,Genbank:NM_099983.2;Name=NM_099983.2;gbkey=mRNA;gene=NAC001;inference=similar to RNA sequence%2C mRNA:INSD:BT001115.1%2CINSD:AF439834.1%2CINSD:AK226863.1;locus_tag=AT1G01010;orig_protein_id=gnl|JCVI|AT1G01010.1;orig_transcript_id=gnl|JCVI|mRNA.AT1G01010.1;product=NAC domain containing protein 1;transcript_id=NM_099983.2
	while (getline(ingff3, line0))
	{
		if (line0.find("#") == 0) //����#��ͷ����
		{
			continue;
		}
		std::vector<std::string> mystr; //�����и���ַ���
		mystr = split(line0, "\t");
		std::cout << mystr[2] << std::endl;
		if (mystr[2] == "gene" || mystr[2] == "Gene" || mystr[2] == "GENE") //��Щggf3�ļ�û��mRNA  �󳦸˾�gff3�ļ���������
		{

			std::vector<std::string> mystr1 = split(mystr[8], ";"); //��ȡ���һ�еķָ��ַ���
			std::vector<std::string> myname = split(mystr1[2], "="); //��ȡname��
			std::vector<std::string> mygeneid = split(mystr1[0], "="); //��ȡ������
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
				trans_info[result1[1]]["Exonarray"].push_back(mystr[4]); //result[0]�Ǳ�ƥ����ַ����������ǲ�����ַ�����1����
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
	ingff3.close(); //��ȡgff3�ļ����
	//����ļ�
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
				if (!(trans_info[iso1->first]["Exonarray"].size()>=1)) //��������������ֵ�Ļ�
				std::cout << "Number " << trans_info[iso1->first]["Cdsarray"][0] << std::endl;
				std::string trans_seq, cds_seq; int trans_len=0, cds_len=0;
				if (!(trans_info[iso1->first]["Exonarray"].size()>=1)) //��������������ֵ�Ļ�
				{
					std::cout<< "WARN: no exon in "<< iso1->first <<std::endl;
					continue;
				}
				std::vector<int> exon=sortstring(trans_info[iso1->first]["Exonarray"]); //�����С����
				std::vector<int> cds =sortstring(trans_info[iso1->first]["Cdsarray"]); //��С��������
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
				//std::vector<std::string>::iterator iter = std::find(trans_info[iso1->first]["strand"].begin(), trans_info[iso1->first]["strand"].end(), '-');//���ص���һ��������ָ��
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
				int N_num = count(trans_seq.begin(), trans_seq.end(), 'N'); //ͳ��������N�ĸ���
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
	time(&start); //��ȡ��ʼ��ʱ��
	using std::endl;
	unsigned long count_del_slash_n = 0;
	unsigned long gene_count = 0; //ͳ�ƻ������
	unsigned long transcript_count = 0; //ͳ��ת¼����
	unsigned long exon_count = 0; //ͳ����������
	unsigned long cds_count = 0;//ͳ��cds
	std::map<std::string, std::string> genomefa; //�洢������fa�ļ�
	std::ifstream ingff3(pathingff3, std::ios::in);  //��gff3�ļ�
	std::ifstream ingenome(pathingenome, std::ios::in); //��genome�ļ�
	//��ʼ��ȡ�������ļ�������map����
	std::string temp;
	while (getline(ingenome,temp,'>'))
	{
		if (temp.length() == 0 || temp.find(0x0A) == 0)   //����������һ���Ƿ��ǿջ��߻س�,���з�����0λ�ã�0x0A=��\n��
		{
			continue; //������ͷ�Ŀ���
		}
		size_t pos_n = temp.find('\n'); //���ҵ�һ�λ��з����ֵ�λ�ã�Ϊ�˵õ������ID;size_t��unsigned int����
		std::string gene_id, gene_seq;
		gene_id.assign(temp, 0, pos_n); //���ǶԵ�,��ȡ��һ������ǰ���ַ�����Ϊgene_id
		if (size_t pos_space = gene_id.find(0x20))  //�ҵ�space �ո�
		{
			gene_id.assign(gene_id, 0, pos_space); //��ȡ
		}
		if (size_t pos_tab = gene_id.find(0x09)) //�ҵ�ˮƽ�Ʊ��
		{
			gene_id.assign(gene_id, 0, pos_tab);
		}
		gene_seq.assign(temp, pos_n + 1, temp.length() - pos_n); //+1��Ϊ��ȥ����β�Ļ��з���  //���Ҳ������(��ȡ����)
		//remove������ͷ�ļ�#include<algorithm>�У����������ͷ�ļ�remove�����Ͳ���ʹ��
		//gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), '\n'), gene_seq.end()); 
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x0A), gene_seq.end()); //ɾ�����з���ͬ��
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x09), gene_seq.end()); //ɾ����ˮƽ�Ʊ����
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x08), gene_seq.end()); //ɾ��backspace
		gene_seq.erase(remove(gene_seq.begin(), gene_seq.end(), 0x20), gene_seq.end()); //ɾ���ո�
		//genomefa[gene_id] = gene_seq;  //fa����������map
		genomefa.insert(std::pair<std::string, std::string>(gene_id, gene_seq)); //����Ҳ��
	}
	printf("fa�����ڴ����\n");
	//fa�����ڴ����

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
		if (line0.find("#") == 0) //����#��ͷ����
		{
			continue;
		}
		std::vector<std::string> mystr;
		mystr = split(line0, "\t");
		if (mystr[2] == "gene" || mystr[2] == "Gene" || mystr[2] == "GENE")
		{
			chrom[mystr[0]]++; //ͳ��Ⱦɫ����Ŀ��һ��Ⱦɫ�����ж��ٸ�����
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
	outlog << "Chromosome numbers:" << chrom.size() << endl; //��ӡȾɫ����Ŀ
	std::map<std::string, int>::iterator ibegin = chrom.begin();
	for (; ibegin != chrom.end(); ibegin++)
	{
		outlog << "Chromosome Name: " << ibegin->first << "\tGene numbers: " << ibegin->second << endl; //��ӡһ��Ⱦɫ�����ж���������
	}
	ingff3.clear(std::ios::goodbit);//��seekg֮ǰ��Ҫ�����������clear�����������ı�������������Ժ�Ϳ�����������seekg�����ˡ�
	ingff3.seekg(std::ios::beg); //begin����д,���ļ�ָ�����»ص�ͷ�� //seekp����������ļ�ָ��λ�ã�seekg�����������ļ�ָ��λ��
	std::string line;
	std::string catcds;
	while (getline(ingff3, line)) //Ĭ�ϻ��з�Ϊһ��
	{
		if (line.find("#") == 0) //����#��ͷ����
		{
			continue;
		}
		std::vector<std::string> mystr;
		mystr = split(line, "\t");
		if (genomefa.count(mystr[0])>0)  //˵���������
		{
			if (mystr[2] == "gene" || mystr[2] == "Gene" || mystr[2] == "GENE")
			{
				std::vector<std::string> mystr1 = split(mystr[8], ";"); //��ȡ���һ�еķָ��ַ���
				std::vector<std::string> myname = split(mystr1[2], "="); //��ȡname��
				std::vector<std::string> mygeneid = split(mystr1[0], "="); //��ȡ������
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
			if (mystr[2] == "CDS" || mystr[2] == "Cds" || mystr[2] == "cds") //��CDSƴ��һ��
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
	time(&end); //��ȡ������ʱ��
	printf("Elapse %d min\n", (unsigned int)(end - start) / 60); //��ӡʱ��
}



