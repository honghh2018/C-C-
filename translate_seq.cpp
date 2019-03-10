#include<iostream>
#include<fstream>
#include<string>
#include<cassert>
#include<stdlib.h>
#include<algorithm> //need reverse function
#include <vector>
#include <stdio.h>
using namespace std; //need declare namespace for string data type 

void readfile(string str); //avow the function
typedef char * pchar;
int main(int argc,pchar argv[])  //array of pointor 
{
	if (argc != 2)
	{
		cout << "your input wrong, please taking overhaul\n";
		exit(0);
	}
	string tst(argv[1]);
	//string tst; //windows pathwar needed double slash
	cout << tst << endl;
	readfile(tst);
	system("pause");
	return 0;
}
/*ATCG
TTGC
CCTG*/
void readfile(string str)
{	
	ifstream infile(str);
	if (!infile.is_open()) 
	{
		cout << "Can not open this file" << endl;
		assert(infile.is_open()); //if did not open file output the wrong information
	}
	string ori;
	cout << "origin string" << endl;
	while (getline(infile, ori))
	{
		cout << ori << endl;
	}
	infile.clear(); //del content
	infile.seekg(0,ios::beg); //file handle of pointor return to start
	cout << "Translate start: " << endl;
	string s;
	while (getline(infile, s)) 
	{
		string::iterator iterA;
		for (iterA = s.begin(); iterA != s.end(); iterA++)
		{
			if (*iterA == 'A' || *iterA=='a')
			{
				*iterA = 'T';
			}else if (*iterA == 'T'|| *iterA == 't')
			{
				*iterA = 'A';
			}else if (*iterA == 'C' || *iterA == 'c')
			{
				*iterA = 'G';
			}else
			{
				*iterA = 'C';
			}
		}
		reverse(s.begin(), s.end());
		cout << "sting: "<<s<< endl;
	
	}
	infile.close();
}

/*Usage:
    through vs 2017 to compile this source code,it since used cmd to run this program
    translate_seq.exe pathway for windows file
*/


