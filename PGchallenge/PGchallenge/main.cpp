#include <windows.h>
#include <stdlib.h>
#include <iostream>
#include <atlstr.h>

#include "timer.h"
#include  "PGAlg.h"


using namespace std;
using namespace TimeCounters;


#define  MAX_TIME  20*60

void main(int iArgCount, char **ppcArgValues)
{
	CTimeCounter  c_time_counter;
	double  d_time_pased;

	
	try
	{
		CPGAlg  c_test;

		c_test.bInitialize("test.txt");

		c_time_counter.vSetStartNow();

		for (int i_iter = 0; d_time_pased < MAX_TIME; i_iter++)
		{
			c_test.vRunIteration();
			cout << "iteration " << i_iter << "current best tree: " << c_test.sGetCurrentBestTree();

			c_time_counter.bGetTimePassed(&d_time_pased);
		}//for (int i_iter = 0; d_time_pased < MAX_TIME; i_iter++)
		


		
	}//try
	catch (exception  *p_ex)
	{
		cout << p_ex->what();
		delete p_ex;
	}//catch (exception  *p_ex)

};//void main()