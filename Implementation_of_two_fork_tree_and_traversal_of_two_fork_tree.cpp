#include<iostream>
#include<string>
#include<stack>

//标准的二叉树存储使链式存储，每一个二叉树的节点是一个结构体
using namespace std;

//中序，前序，后序，递归遍历和非递归遍历的实现，非递归遍历需要借助栈和循环
//二叉树的遍历每个节点只访问一次。
//时间效率是O（n）
//空间效率也是O（n）
//树深为k的递归遍历需要k+1个辅助单元

//增强版的遍历，中序加前序遍历的算法。

//中序是左右中，先序是中左右，后序是左右中。

//



struct Mystruct
{   //这里的结构体变量直接定义就初始化，省去了构造函数中的初始化（非常方便）
	int Nodedata=0;  //这个代表二叉树的节点（数据域用来存储数据的）
	//一个二叉树有左节点和右节点
	//C++结构体中的成员指针能给默认值初始化nullptr,防止野指针
	Mystruct *pleft=nullptr; //一个结构体指针，代表左节点（指针域）
	Mystruct *pright=nullptr; //代表右节点（指针域）
}BTrue,*pBTrue;  //BTrue是一个实例化的对象结构体，pBTrue是一个结构体指针

//写一个函数实现打印二叉树容器的数据到屏幕上
void show(Mystruct *pRoot,int n)  //传递外部的二叉树节点指针（不是是根节点指针）
{
	//使用递归呢？递归比栈好写
	if (pRoot == nullptr)  //如果等于空，就直接返回
	{
		return;   //直接返回一个空，什么都不动作的意思
	}
	else
	{
		//先递归左边，再右边（对于一个三个节点的树），n是用来打印空格的cout << " "; 
		show(pRoot->pleft, n + 1);  //这里的pRoot->pright也不影响遍历二叉树
		for (int i = 0; i < n; i++)  //每一行打印一排空格
		{
			cout << "  "; //打印一个空格
		}

		//打印数据
		cout << pRoot->Nodedata << endl; //打印数据

		//上面递归左边，这里递归右边
		show(pRoot->pright, n + 1);
	}

}


//写一个中 序遍历来遍历二叉树，对于中序遍历只需要一个节点就行，函数的形式参数一个就行

void zhongxubianli(Mystruct *proot)  //函数的副本机制导致的proot每次递归等于proot->pleft，然后proot->pleft又等于proot->pleft,不断的递归调用，导致不断的循环至二叉树的叶子节点
{   //中间节点（先从根节点开始，只要不等于nullptr,就一直往纵深递归，一直递归到叶子节点，它等于nullptr后返回）
	if (proot != nullptr)   //等于nullptr就不需要动作了，不等于就需要动作
	{  //中序遍历二叉树的访问顺序是，左边，中间，右边
		if (proot->pleft != nullptr)  //不等于空就继续
		{  //左节点
			zhongxubianli(proot->pleft);//proot->pleft一直往下纵深递归，知道if中的proot->pleft != nullptr判断为空跳出循环，下同
		}  //cout << " " << proot->Nodedata << endl;前序遍历和后序遍历，只要把zhongxubianl中的if语句调用的顺序改一下就可以了（放到前面就是前序遍历，放到后面就后序遍历）
		//中间把数据打印出来
		cout << " " << proot->Nodedata << endl;
		if (proot->pright != nullptr)
		{  //右边节点
			zhongxubianli(proot->pright);
		}
	}


}









//整个二叉树的操作基本是递归完成的
//非递归模式搞定二叉树的遍历（难度比较大）
//使用栈完成遍历二叉树,中序遍历
void stackzhong(Mystruct *proot)
{
	//构建一个栈
	//stack<Mystruct> mystack; //创建我们自己定义的结构体的数据类型的栈
	Mystruct *pcurrent = proot;//记录根节点
	Mystruct *mystack[100];//预留一百个节点(创建一个指针数组)
	int top = 0; //栈顶元素为0下标，记录一下下标，和数组有点区别
	while (top!=0||pcurrent!=nullptr)   //栈顶元素下标为0，说明为空栈，并且节点不能为空
	{
		while (pcurrent != nullptr)  //如果节点不为空，那么我们就压栈
		{
			mystack[top++] = pcurrent;   //top++是先引用再加加
			pcurrent = pcurrent->pleft;//因为中序遍历，我们先左边，不断的纵深递归
		}

		if (top > 0)
		{
			top--;//数据出栈
			pcurrent = mystack[top];  //栈顶元素
			cout << " " << pcurrent->Nodedata << endl;
			//访问右边节点
			pcurrent = pcurrent->pright;
		}

	}


}

//优化一下stackzhong,如下：stackzhongA

void stackzhongA(Mystruct *proot)
{
	Mystruct *pcurrent = proot;
	stack<Mystruct*> mystack;  //创建一个栈对象
	while (!mystack.empty() || pcurrent != nullptr)   //栈非空!mystack.empty() 就继续
	{
		while (pcurrent != nullptr)  //如果节点不为空，那么我们就压栈
		{
			mystack.push(pcurrent);   //将节点指针进栈
			pcurrent = pcurrent->pleft;//因为中序遍历，我们先左边，不断的纵深递归
			//将左边的节点全部入栈
		}

		if (!mystack.empty())  //非空的情况下，数据出栈
		{
			pcurrent = mystack.top();  //不断的循环，因为栈顶元素不断的pop出栈，这里就循环遍历了栈中的元素，然后求出它的节点的数据域
			cout << " " << pcurrent->Nodedata << endl;  //数据出栈
			mystack.pop();//弹出栈的元素
			
			//访问右边节点
			pcurrent = pcurrent->pright;
		}

	}


}

// 写一个算法实现求出树的叶子。
//用递归求叶子节点个数
int getleaf_num(Mystruct *proot)  //传入的形式参数是mystruct结构体类型的指针
{

	if (proot == nullptr)  //意外情况
	{
		return 0; //返回0
	}
	if (proot->pleft == nullptr && proot->pright == nullptr)  //这个是终止条件，也就是说这个树只有一个根节点，没有左右节点（子叶）
	{  //左右节点为空，就返回1个叶子节点
		return 1; //返回一个节点（有多少个函数返回值，就进行返回值的相加）
	}
	//它会一直递归下去，因为函数形式参数的副本机制，一直递归到叶子，后再一层一层的返回
	int left = getleaf_num(proot->pleft);  //左叶子
	int right = getleaf_num(proot->pright);//右叶子
	//返回
	return left + right;
}

//求一个完全二叉树的深度，这里不是满二叉树是完全二叉树，如下的s1节点的二叉树
int getheight(Mystruct *proot)
{
	int height = 0;  //深度
	int left = 0;  //左
	int right = 0;//右
	if (proot == nullptr)  //如果节点为空(终止条件我们返回0)
	{
		return 0; //或者return 0或者return height
	}
	//实现函数自己调用自己就是递归
	left = getheight(proot->pleft);
	right = getheight(proot->pright);

	//深度取决于左边或者右边最大一个节点数
	//返回left和right后取它们的最大值，使用三目运算符
	height = left > right ? left : right; //三目运算符，如果left大于right那么返回left，否则返回right
	return right + 1; //
}

//递归总结
//递归其实不难就是靠递归来实现操作，递归一个函数，等它碰到条件被破坏后，结束递归，不断的返回，返回中下面的语句就是递归函数实现操作的关键语句了。或者在递归函数上方，在递归的时候不断的执行这个语句完成操作

//判断一个二叉树是完全二叉树还是满二叉树，就是遍历每一个节点

//按照队列的方式来访问所有节点，遍历所有节点
//按层次遍历
void ceng(Mystruct *proot)
{  //这里的进队是12345678，出队也是12345678
	if (proot == nullptr)  //如果根节点为空，直接返回，什么都不动作
	{
		return; //下面代码不执行了，提前结束函数
	}
	Mystruct *myq[100]; //指针数组
	//队列的头部和尾部，如下：
	int tou = 0;
	int wei = 0;
	//一开始初始化为空
	Mystruct *pcurrent = nullptr;  //当前位置指针
	//wei++;  //尾部递增
	//如果头部和尾部重合说明队列没有元素，|头   |尾，那么||这里的头部和尾部重合了，说明队列没有元素，这里|   |头部和尾部没有重合说明队列里面有元素，尾部往后移动，头部不动是入队，头部往后动了是出队，尾部往后动了是入队
	myq[wei++] = proot;//（根节点1入队，下面的if判断是左节点和右节点进队）  //直接下标wei++也行，存入队列的第一个二叉树节点
	while (tou != wei)  //tou不等于wei就继续，头部不等于尾部就继续
	{ //tou++是出队和栈的pop函数一样，wei++是入队，和栈的push函数一样
		pcurrent = myq[tou];//当前位置的指针等于头部
		tou++; //1出队//头部递增////如果头部和尾部重合说明队列没有元素，|头   |尾，那么||这里的头部和尾部重合了，说明队列没有元素，这里|   |头部和尾部没有重合说明队列里面有元素，尾部往后移动，头部不动是入队，头部往后动了是出队，尾部往后动了是入队
		//打印当前指针pcurrent指向的数据
		cout << pcurrent->Nodedata << endl;
		//if (pcurrent->pleft)  //判断左节点不为空就继续
		if (pcurrent->pleft != nullptr)  //这个和上面的写法意义是一样的，不为空就继续
		{
			//需要入队，左入队
			myq[wei++] = pcurrent->pleft; //指针数组元素等于给定地址

		}
		if (pcurrent->pright)  //右节点不等于nullptr就继续
		{
			//右入队
			myq[wei++] = pcurrent->pright; //右节点入队
		}
		
	}
}

//二叉树如何找父节点和根节点？
//寻找父节点函数
int get_father(Mystruct *pRoot, int num)
{  

	if (pRoot == nullptr)
	{
		return 0; //如果pRoot等于空就返回零
	}
	//非空就返回父节点的数据，这是左节点
	if (pRoot->pleft != nullptr&&pRoot->pleft->Nodedata == num)
	{
		return pRoot->Nodedata;  
	}
	//这是右节点
	if (pRoot->pright != nullptr&&pRoot->pright->Nodedata == num)  ////寻找二叉树的父节点函数调用,这里是寻找num的父节点，num是二叉树的数据元素，值为num不是二叉树的下标
	{
		return pRoot->Nodedata;
	}
	//判断父节点需要递归
	get_father(pRoot->pleft, num);
	get_father(pRoot->pright, num); //不断的向纵深递归到叶子节点的左右节点
}

//如果知道右兄弟，我们要找它的左兄弟如何设计算法呢？

int get_left(Mystruct *pRoot, int num)
{
	if (pRoot == nullptr)
	{
		return 0; //如果pRoot等于空就返回零
	}
	//非空就返回父节点的数据，这是右节点
	if (pRoot->pright != nullptr&&pRoot->pright->Nodedata == num)  //知道右节点，寻找左节点
	{
		if (pRoot->pleft!=nullptr)  //在左节点不为空的情况下
		return pRoot->pleft->Nodedata;  //返回左节点数据，求左兄弟
	}
	
	//知道右节点寻找左节点需要递归
	get_left(pRoot->pleft, num);
	get_left(pRoot->pright, num);


}

//知道一个二叉树的左节点求右节点数据
int get_right(Mystruct *pRoot, int num)
{
	if (pRoot == nullptr)
	{
		return 0;// 直接返回什么都不动作
	}
	if (pRoot->pright != nullptr&&pRoot->pleft->Nodedata == num)
	{
		return pRoot->pright->Nodedata;
	}

	//实现递归
	get_right(pRoot->pleft, num);
	get_right(pRoot->pright, num);
}

//实现二叉树的插入，删除，修改，查找，排序
//难点是插入和删除和排序
//二叉树需要平衡，不能只插左节点或者只插右节点，只能是用平衡二叉树
//左右插才能平衡，不能只插一边。在进行插入算法，并且是平衡的不能只插入左边或者右边，这就非常复杂了（插入复杂度本身不高是O（1））
//插入和删除依赖于查找

//红黑树就是平衡二叉树。

int findmax(Mystruct *proot)  //找出二叉树最大的元素
{
	int max = -99999; //一开始我们初始化它为一个比较小的值
   //我们使用中序遍历，我们使用的是结合栈的
	//构建一个栈
	//stack<Mystruct> mystack; //创建我们自己定义的结构体的数据类型的栈
	Mystruct *pcurrent = proot;//记录根节点
	Mystruct *mystack[100];//预留一百个节点(创建一个指针数组)
	int top = 0; //栈顶元素为0下标，记录一下下标，和数组有点区别
	while (top != 0 || pcurrent != nullptr)   //栈顶元素下标为0，说明为空栈，并且节点不能为空
	{
		while (pcurrent != nullptr)  //如果节点不为空，那么我们就压栈
		{  //数据进栈的过程
			mystack[top++] = pcurrent;   //top++是先引用再加加
			pcurrent = pcurrent->pleft;//因为中序遍历，我们先左边，不断的纵深递归
		}

		if (top > 0)
		{
			top--;//数据出栈
			pcurrent = mystack[top];  //栈顶元素
			//cout << " " << pcurrent->Nodedata << endl;
			if (max < pcurrent->Nodedata)
			{
				max = pcurrent->Nodedata; //如果max小于pcurrent->Nodedata那么用max接收Nodedata的值
			}
			//访问右边节点
			pcurrent = pcurrent->pright;
		}

	}

	return max; //返回最大值
}


//实现二叉树的元素插入
//平衡数里面删除一个数据非常的麻烦，因为删除一个元素后还要保存平衡

//这里我们实现插入，插入要有一个返回值
Mystruct *insertnode(Mystruct *proot, int num)  //插入一个节点
{  //如果这个二叉树的根节点为空，我们直接开辟堆内存和插入元素数据
	if (proot == nullptr)  //如果头指针为空
	{
		Mystruct *pnew = new Mystruct;  //开辟堆内存空间
		pnew->Nodedata = num; //初始化
		proot = pnew;
	}
	else if (num<=proot->Nodedata)  //为了实现插入二叉树的数据按顺序排序，这里需要判断
	{   //小于等于放左边
		proot->pleft = insertnode(proot->pleft, num); //递归实现
	}
	else
	{  //大于等于放右边
		proot->pright = insertnode(proot->pright, num); //递归实现
	}
	return proot;
}

//删除叶子节点很简单，但是删除中间节点非常苦难，还需要将后面和前面的数据连接起来
//又要删除二叉树的中间节点，又要移动数据，又要平衡二叉树，非常难的
//一般面试不会考你二叉树的平衡，但是会考你两个二叉树的祖先

//如何判定一个二叉树是完全二叉树吗？

//线索二叉树，一般二叉树只能访问左右节点



//树和森林（小公司不需要你现场写二叉树，大公司需要）
//平衡树，不用set的情况，非常复杂。
//有一个清华小伙总结的数据结构笔记是一个word
//霍夫曼树是一个压缩的树存储。




//实现插入的main
void main()
{
	//指针在C++中必须初始化
	Mystruct *pRoot = nullptr; //创建一个指向二叉树的根节点的指针pRoot
	for (int i = 6; i < 10; i++) //左插5个
	{
		pRoot = insertnode(pRoot, i); //插入1到10的数据元素
	}
	//实现右插（这里是为了平衡二叉树，左插5个，右插5个），但是还是不平衡，类似于链表
	for (int i = 5; i > 0; i--) //右插5个
	{
		pRoot = insertnode(pRoot, i); //插入1到10的数据元素
	}


	//显示数据
	show(pRoot,1);  //1是空格层次



	cin.get();
}



void main11111()
{
	//构建一个简单的二叉树
	//创建一个二叉树的根节点，使用指针创建
	Mystruct *pRoot; //代表二叉树的根节点,pRoot我们让它等于s1的地址&s1
	//这是一个指针指向根节点。
	//可以根据数组初始化，这个二叉树

	//实例化多个结构体节点对象，如下：
	Mystruct s1;  //每个对象是一个二叉树的节点
	Mystruct s2;
	Mystruct s3;
	Mystruct s4;
	Mystruct s5;
	Mystruct s6;
	Mystruct s7;
	//
	Mystruct s8;
	//代表二叉树的根节点,pRoot我们让它等于s1的地址&s1
	pRoot = &s1;


	//初始化数据，并将这些结构体节点连接起来变成二叉树
	//初始化结构体数据
	s1.Nodedata = 1;  //相当于链表节点的数据域
	s1.pleft = &s2; //s1的左指针等于s2结构体节点，这就连接起来了
	s1.pright = &s3; //取右节点结构体地址
	s2.pleft = &s4;
	s2.pright = &s5;
	s3.pleft = &s6;
	s3.pright = &s7;
	s4.pleft = &s8;  //s4节点的pleft等于&s8
	//初始化数据域
	s2.Nodedata = 2;
	s3.Nodedata = 3;
	s4.Nodedata = 4;
	s5.Nodedata = 5;
	s6.Nodedata = 16; //最大值findmax函数
	s7.Nodedata = 7;
	//
	s8.Nodedata = 8;

	//能创建数组来实现二叉树数据进行初始化（但是比较麻烦）


	//显示二叉树容器的数据，使用递归函数
	//show(pRoot, 1);  //这里传入1是打印1个空格
	//调用中序调用函数
	//zhongxubianli(pRoot);
	//中序遍历能实现二叉树的查找

	//使用栈实现中序遍历二叉树
	//stackzhong(pRoot);
	//用优化的栈实现遍历二叉树，中序遍历
	//stackzhongA(pRoot);

	//求二叉树的叶子节点个数,这里结果是4个叶子
	//cout<<getleaf_num(pRoot)<<" ";
	//打印一个二叉树的深度
	//cout << getheight(pRoot) << endl;
	
	//层次感遍历
	//ceng(pRoot);  //结果是12345678
	//

	//寻找二叉树的父节点函数调用,这里是寻找3的父节点，3是二叉树的数据元素，值为3不是二叉树的下标
	//cout << get_father(pRoot, 3) << endl;
	//寻找二叉树数据元素为7的父节点，必须是满二叉树
	//cout << get_father(pRoot, 7) << endl;

	//知道右兄弟，寻找左兄弟,知道一个右兄弟3寻找它的左兄弟，如果是知道一个2，2已经是左兄弟了，它没有左兄弟
	//cout << get_left(pRoot,3) << endl;
	
	//知道左兄弟求右兄弟,知道左兄弟是2，求右兄弟
	//cout << get_right(pRoot, 2) << endl;

	//cout << get_right(pRoot, 6) << endl;

	//找出二叉树的最大值,并打印出来，这里最大是8
	cout<<findmax(pRoot) << endl;


	cin.get();
}














///我们用难度很高的初始化二叉树的方法，我们使用数组初始化它
void main2222()
{
	//定义一个mystruct类型的数组
	Mystruct struct_array[100];
	//首先定义根节点
	Mystruct *pRoot;
	pRoot = struct_array;//这里是第一个元素的首地址，数组名不用取地址符号就是第一个元素的首地址
	
	//给这个mystruct类型数组赋值,for循环，这个mystruct数组每一个元素是一个结构体
	for (int i = 0; i < 100; i++)
	{
		struct_array[i].Nodedata = i; //对结构体数组中的结构体成员赋值
	}
	
	//推理一下如果第一个节点的数据域是0，那么它存储数据域为1和2的两个节点的首地址（指针域）
	//数据域为1的节点存储数据域为3和4节点的首地址
	//那么它们的关系是，从零开始，节点数为n的那么它下一个节点数是2*n+1和2*n+2,其实也不叫数据域了，就是一个和数组一样的二叉树的下标，（二叉树的下标也是从零开始）
	//那么，推理一下如果n=1,那么n下面的两个几点的下标是2*n+1和2*n+2分别为2*1+1=3和2*1+2=4

	//这时候用这个结构体数组初始化二叉树，struct类型pRoot
	//2*i+2<99要满足这个条件，避免数组越界
	//for (int i = 0; i <= (99-2) / 2;i++)  //因为左右节点都是2*n+1和2*n+2,n是一个二叉树的下标，和数组下标一样，因此这里的100必须除以2，才能表示二叉树的下标，否则2*100越界了。2*100/2正好是100个数组下标初始化100个二叉树节点下标
	for (int i = 0; i <= 50; i++)  //这里显示的是从0-99才正确
	{    //0-98就是99个数，少一个数；99-1是越界了，99-2少两个数，如何让它一个数都不少呢？使用if判断，并把for的循环条件改成i<50
		if (i <= (99 - 1) / 2)  //98除以2
		{
			//先初始化左节点
			struct_array[i].pleft = &struct_array[2 * i + 1];
		}
		if (i <= (99 - 2) / 2)  //97除以2
		{
			//再初始化右节点（其实中间节点也是上一个节点的左节点或者右节点，那么所有的节点都能看成是左右节点而没有中节点）
			struct_array[i].pright = &struct_array[2 * i + 2]; //取数组中的每一个元素的地址给二叉树指针域赋值
		}
		
		
	}
	//显示这个二叉树，调用show函数
	show(pRoot,1);


	cin.get();
}
