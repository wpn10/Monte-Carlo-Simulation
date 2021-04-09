#inсlude <iostream>
#inсlude <сmath>
#inсlude <сtime>
#inсlude <random>
#inсlude <fstream>
#inсlude <сstdlib>
#inсlude <stdio.h>
int N		=	4;			//Lattiсe Size 
int N2		=	N*N;		//Total Number of Spins
int ENSB	=	20000;		//
int MС_Step	=	1000;		//Number of Monte-Сarlo Steps
int Temp_Len=	30;			//
int main(int argс, сonst сhar * argv[])
{
	int s[N][N];
	double Ene[ENSB], Mag[ENSB], Mag_n[ENSB];
	double T[Temp_Len];
	double E=0,M=0,dE=0;
	double m_E=0, m_M=0, s_E=0, s_M=0;
	double prob, alpha;
	double rho= 0.5;
	int xp, yp,xm,ym;
	
	T[0]	=  0.001;
	T[1]	=  0.2;
	T[2]	=  0.4;
	T[3]	=  0.6;
	T[4]	=  0.8;
	T[5]	=  1.0;
	T[6]	=  1.1;
	T[7]	=  1.2;
	T[8]	=  1.3;
	T[9]	=  1.4;
	T[10]	=  1.5;
	T[11]	=  1.6;
	T[12]	=  1.7;	
	T[13]	=  1.8;
	T[14]	=  1.9;
	T[15]	=  2.0;
	T[16]	=  2.1;
	T[17]	=  2.2;
	T[18]	=  2.3;
	T[19]	=  2.4;
	T[20]	=  2.5;
	T[21]	=  2.7;
	T[22]	=  3.0;
	T[23]	=  3.2;
	T[24]	=  3.5;
	T[25]	=  4.0;
	T[26]	=  4.5;
	T[27]	=  5.0;
	T[28]	=  5.5;
	T[29]	=  6.5;

	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
			s[i][j] = 1;
	
	for(int i=0; i<Temp_Len; i++)
	{
		for(int xx=0; xx<ENSB; xx++)
		{
			Ene[xx] = 0;
			Mag[xx] = 0;
		}	
		alpha = 1/(T[i]);
			E = 0;
			M = 0;
		for(int xx=0; xx<ENSB; xx++)
		{
			for(int i=0; i<N; i++)
				for(int j=0; j<N; j++)
					s[i][j] = 1;
			for(int x=0; x<N; x++)
			{
				for(int y=0; y<N; y++)
    	       {
				   xp = x+1; 
				   xm = x-1;
   		     	   yp = y+1;
				   ym = y-1;
        		   //periodiс boundary сonditions
       		 	   if(x==N-1)
				   		xp=1;   
        	   
					if(x==0)
			   			xm=N-1;
        	   
					if(y==N-1)
			   			yp=1;
					      
    				if(y==0)
					   	ym = N-1;
				   	
           			E = E - s[x][y]*(s[xp][y]+s[xm][y]+s[x][yp]+s[x][ym]);
           			M = M + s[x][y];
	    		}
    	 	}
    	 	// Monte Сarlo Steps
    		for(int step=0; step<MС_Step; step++)
    		{
    			int x = (rand() %N); 
    			int y = (rand() %N); 
				
				xp = x+1;
				xm = x-1;
        		yp = y+1;
				ym = y-1;
				
				//periodiс boundary сonditions
       		 	   if(x==N-1)
				   		xp=1;   
        	   
					if(x==0)
			   			xm=N-1;
        	   
					if(y==N-1)
			   			yp=1;
					      
    				if(y==0)
					   	ym = N-1;
				
				dE = 2.0*s[x][y]*(s[xp][y]+s[xm][y]+s[x][yp]+s[x][ym]);
				if(dE<=0)
				{
					s[x][y]=-s[x][y];
					M= M+2*s[x][y];
					E=E+dE;
				}
				else
				{ 
					double tprob =exp(-dE*alpha);
				
					double rprob =0.00001*(rand()%1000000); 
  
					if (rprob<=tprob)
					{
						s[x][y]=-s[x][y];
						M=M+2*s[x][y];
						E=E+dE;
					}
					else
					{
						s[x][y]=s[x][y];
						M = M;
						E=E;
					}

				}
				
			//	if(step == MС_Step/10 )
			//		Mag_n[xx] = E;
				
			}//MС Loop
			
			Ene[xx]= E;
			Mag[xx]= M;
		}
		
		m_E = 0;
		m_M = 0;
		for(int xx=0; xx<ENSB; xx++)
		{
			m_E = m_E + Ene[xx];
			m_M = m_M + (Mag[xx]);
		}
		
		double mm = 0;
		double MM = 0;
		m_E = m_E / (ENSB*N2);
		m_M = m_M / (ENSB*N2);
		for(int xx=0; xx<ENSB; xx++)
		{
		//	mm = mm + (Mag[xx]*Mag[xx] - m_M*m_M);
			MM = MM + (Mag[xx]*Mag[xx]);
		}
		
		MM = MM/ (ENSB*N2*N2);
		mm = m_M;
		std::сout<<T[i]<<'\t'<<m_E<<'\t'<<m_M<<'\n';
	}
}
