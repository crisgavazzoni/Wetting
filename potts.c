/*********************************************************************************************
*																							 *
*          Programa gota (oleo e agua) sobre uma superfície de pilares inclinada             *
*			                Cristina Gavazoni, Outubro 2024									 *
*																							 *
*********************************************************************************************/
/* Não está implementado a continuação de runs antigos!!!*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mc.h"
#include "pointers.h"
#ifdef USEGFX
#include "board3d.h"
#endif

/*********************************************************************************************
*                                 Declaração das funções                                     *
*********************************************************************************************/

void create_pointers(void);
void free_pointers(void);
void openfiles(void);
void initialization(void);
void initial_state(int *s, int l, int rg, int initialstate);
void xyz(int *s, int l, int m);
void save_conf(int num_steps, int iout);
void dynamics (int *s,int num_steps);
void flip_spin (long site, long new_s);
void add_to_interface(int site);
void remove_from_interface(int site);
void measure_angle(int num_steps);
double calculate_energy(int num_steps);

/*********************************************************************************************
*                       Declarando parâmetros da simulação - técnicos                        *
*********************************************************************************************/

#define mc_steps   	       10000  // Número de passos de MC totais
#define n_mesure       	   1000   // Intervalo para salvar medidas

#define temp           	   13.0  // Temperatura
#define kB             	   1.0  // Constante de Boltzman

#define t_neigh        	   26   // Número de vizinhos
#define t_close_neigh  	   18   // Número de vizinhos próximos 

#define NUM_CONF            1 // DEIXAR IGUAL A 1 !!

/*********************************************************************************************
*                    Declarando Parâmetros constantes durante a simulação                    *
*********************************************************************************************/

#define frac_h 				1.0 // Se > 1 a gota flutua sobre a superfície
#define h_base       		3 //Altura da Base

#define GW                  0.0001 // Massa agua em um sítio
#define GO                  0.000077 // Mass óleo em um sítio
#define Lambda_w            0.01 // Parâmetro de volume
#define Lambda_o            0.01 // Parâmetro de volume

#define eps_SG              0.96 // E superficial S-G

#define eps_WG              2.70 // E superficial L-G agua
#define cos_theta_w         (cos(111.0*M_PI/180.0))
#define eps_SW              (eps_SG-eps_WG*cos_theta_w)// E superficial S-L agua


#define eps_OG              1.04 // E superficial L-G oleo
#define cos_theta_O         (cos(53.0*M_PI/180.0))
#define eps_SO              (eps_SG-eps_OG*cos_theta_O)// E superficial S-L oleo

#define eps_WO              2.06 // E superficial óleo-água

/*********************************************************************************************
*                                   Declarando das variáveis                                 *
*********************************************************************************************/

unsigned long identifier;

//--------------------------------------------------------------------------------------------
// Variáveis da configuração inicial da gota
//--------------------------------------------------------------------------------------------

int l, l2, l3, lm;  			// Tamanho do sistema e multiplos
int rg, rg2, rg3;   			// Tamanho da gota e multiplos
int h, w, a, d;     			// Altura, largura, separaco pilar e d=w+a
int initialstate;  				// Estado inicial: W ou CB
double Omg;						// angulo de inclinacao

double v0, v0_w, v0_o; 			// volumes inicial
int t_vol, t_vol_w, t_vol_o; 	// volumes desejados
double f_w, f_o;                // Frações de água e óleo 

double T;

char CI[4]; 					// Guarda o nome do arquivo

//--------------------------------------------------------------------------------------------
// Ponteiros
//--------------------------------------------------------------------------------------------

int *s ; 										// Spins
int *nv_cw,*nv_co, *nv_w, *nv_o, *nv_s, *nv_g; 	// Números de vizinhos de cada tipo
int *inter_pos, *w_inter, *inter_mtx; 			// Sítios na interface

//--------------------------------------------------------------------------------------------
// Variaveis usadas na simulação
//--------------------------------------------------------------------------------------------

int n_w, n_o, n_s; 				//número de sítios água e óleo
int int_label; 					//label do sítio na interface
int vol, vol_w, vol_o; 			// volumes calculados
double Gw, Go;					// Termos gravitacionais

//--------------------------------------------------------------------------------------------
// Nomes dos arquivos de entrada e saida 
//--------------------------------------------------------------------------------------------

char output_file1[100],input_file[100]; 						//Nome dos arquivos
FILE *fconf,*fp1,*finit,*fxyz,*flast,*fp_in,*fout, *fener, *fsol; // Localizadores dos arquivos

//--------------------------------------------------------------------------------------------
// Estatística
//--------------------------------------------------------------------------------------------

long int trial,accepted,acc_interface,acc_teste;

//--------------------------------------------------------------------------------------------
// Coisas para multiplas simulações
//--------------------------------------------------------------------------------------------

int cont; //controla o numero de continuacoes de uma dada simualcao
int file_in=0,confstep,iconf;
long int tempo_in=0;
unsigned int numsteps;
char lixo;

//--------------------------------------------------------------------------------------------
// Variaveis dos observáveis
//--------------------------------------------------------------------------------------------

int drop_height,drop_bottom;
int center_x,center_y;
double Brx,Bry, Bax, Bay;
double Hrx,Hry, Hax, Hay;
double Rrx,Rry, Rax, Ray;
double Trx,Try, Tax, Tay;

/*********************************************************************************************
==============================================================================================
=                                      PROGRAMA PRINCIPAL                                    =
==============================================================================================
*********************************************************************************************/
int main(int argc,char *argv[])
{

	int i;

//--------------------------------------------------------------------------------------------
//  Lendo os parâmetros de entrada e definindo parametros da conf inicial
//--------------------------------------------------------------------------------------------

	identifier = 0;

	if(argc==19) 
	{
   		 for (i=1;i<argc;i++) 
		{
      		if (!strcmp(argv[i],"-L"))        l=atoi(argv[++i]);
      		else if (!strcmp(argv[i],"-R"))   rg=atoi(argv[++i]);
      		else if (!strcmp(argv[i],"-Omg")) Omg=atof(argv[++i]);
      		else if (!strcmp(argv[i],"-fo"))  f_o=atof(argv[++i]);
      		else if (!strcmp(argv[i],"-CI"))  initialstate=atoi(argv[++i]);
      		else if (!strcmp(argv[i],"-w"))   w=atoi(argv[++i]);
      		else if (!strcmp(argv[i],"-h"))   h=atoi(argv[++i]);
      		else if (!strcmp(argv[i],"-a"))   a=atoi(argv[++i]);
      		else if (!strcmp(argv[i],"-s"))   seed=atol(argv[++i]);
      		else {
				fprintf(stderr,"Error.  Argument '%s' is not recognized.\n",argv[i]);
				exit(-1);
      			}// fim do ELSE
    	} // fim do FOR

    	if(initialstate==1) {
    	  sprintf(CI,"CB_");
    	}
    	else if(initialstate==2) {
    	  sprintf(CI,"WE_");
    	}
    	else {
		exit(-1);
    	}

		sprintf(output_file1,"%sgota_3d_L_%d_R_%d_Omg_%.1f_w_%d_h_%d_a_%d.out",CI,l,rg, Omg,w,h,a);	
		fout = fopen(output_file1,"w");	

		fflush(fout);

		d=a+w;

		rg2 = rg*rg;
		rg3 = rg2*rg;
	
		l2 = l*l;
		l3 = l2*l;
		lm = l/2; 

		f_w =  1.0 - f_o;

		v0 = (4.0*M_PI*rg3/3.0);
		v0_w = f_w*v0;
		v0_o = f_o*v0;

		t_vol = (int) v0;
		t_vol_w = (int) v0_w;
		t_vol_o = (int) v0_o;

		T = temp;

		fprintf(fout,"=====================================================================\n");
		fprintf(fout,"=                 Parâmetros Iniciais do Sistema                    =\n");
		fprintf(fout,"=====================================================================\n");

		fprintf(fout,"L      : %d\n",l);
		fprintf(fout,"R      : %d\n",rg);
    	fprintf(fout,"Omg    : %.1f\n", Omg);
		fprintf(fout,"V      : %d\n",t_vol);
    	fprintf(fout,"V_w    : %d\n",t_vol_w);
    	fprintf(fout,"V_o    : %d\n",t_vol_o);
    	fprintf(fout,"fw     : %.2f\n",f_w);
    	fprintf(fout,"fo     : %.2f\n",f_o);
    	fprintf(fout,"h      : %d\n",h);
    	fprintf(fout,"w      : %d\n",w);
    	fprintf(fout,"a      : %d\n",a);
    	fprintf(fout,"eps_SG : %.2f\n", eps_SG);
    	fprintf(fout,"eps_WO : %.2f\n", eps_WO);
    	fprintf(fout,"eps_SW : %.2f\n", eps_SW);
    	fprintf(fout,"eps_WG : %.2f\n", eps_WG);
    	fprintf(fout,"eps_SO : %.2f\n", eps_SO);
    	fprintf(fout,"eps_OG : %.2f\n", eps_OG);	
    	fprintf(fout,"CI     : %s\n",CI);
    	fprintf(fout,"Seed   : %ld\n",seed);
		fprintf(fout,"=====================================================================\n");
		fprintf(fout,"\n");

	} // fim IF 

	else {
    	fprintf(stderr,"Error.     Invalid outout. Use:\n");
    	exit(-1);
  	} //Fim o ELSE

//--------------------------------------------------------------------------------------------
//  Escolhendo uma semente aleatóriamente
//--------------------------------------------------------------------------------------------

	identifier = start_randomic_seed_ext(seed);
	if(identifier!=seed) {printf("id: %ld s: %ld\n",identifier,seed);exit(-1);}

//--------------------------------------------------------------------------------------------
//  Criando Ponteiros
//--------------------------------------------------------------------------------------------

	create_pointers();

//--------------------------------------------------------------------------------------------
//  Abrindo arquivos
//--------------------------------------------------------------------------------------------

	openfiles();

//--------------------------------------------------------------------------------------------
//  Gerando a Configuração inicial
//--------------------------------------------------------------------------------------------

	initialization();
	xyz(s,l,0);
	save_conf(0,0);

//--------------------------------------------------------------------------------------------
//  Estatística
//--------------------------------------------------------------------------------------------

	trial=0; accepted=0;acc_interface=0;acc_teste=0;

//--------------------------------------------------------------------------------------------
//  Dinâmica e medida dos observáveis
//--------------------------------------------------------------------------------------------


	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"=                       Energias do Sistema                         =\n");
	fprintf(fout,"=====================================================================\n");

	for (i=0; i<=mc_steps; ++i) 
	{
		Gw=GW;
		Go=GO;

		dynamics(s,i);

		if(i%n_mesure==0) 
		{
			measure_angle(i);
			save_conf(i,0);
			
		}//fim do IF

	}//Fim do FOR

	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"\n");

	save_conf(mc_steps,1);
	xyz(s,l,1);

//--------------------------------------------------------------------------------------------
//  Estatística
//--------------------------------------------------------------------------------------------

	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"=                         Fim da simulação                          =\n");
	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"Número de partículas de água    : %d\n",vol_w);
	fprintf(fout,"Número de partículas de óleo    : %d\n",vol_o);
	fprintf(fout,"Número de partículas de líquido : %d\n",vol);
	fprintf(fout,"Número de sítios na interface   : %d\n",int_label);
	fprintf(fout,"Fração de óleo                  : %f\n",f_o);
	fprintf(fout,"Fração de água                  : %f\n",f_w);
	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"\n");
	fflush(stdout);
	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"=                            Estatísticas                           =\n");
	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"Trial: %ld\n",trial);
	fprintf(fout,"Accepted: %ld     Ratio: %f\n",accepted,(float)accepted/trial);
	fprintf(fout,"Interface: %ld    Ratio: %f\n",acc_interface,(float)acc_interface/trial);
	fprintf(fout,"Teste: %ld        Ratio: %f\n",acc_teste,(float)acc_teste/trial);
	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"\n\n\n");


//--------------------------------------------------------------------------------------------
//  Encerrando
//--------------------------------------------------------------------------------------------

 	free_pointers();
 	fclose(fout);
 	fclose(fp1);

} //Encerra o Programa
/********************************************************************************************/


/*********************************************************************************************
==============================================================================================
=                                    Subrotinas e Funções                                    =
==============================================================================================
*********************************************************************************************/

/*********************************************************************************************
*                               Rotina que abre os ponteiros                                 *
*********************************************************************************************/
void create_pointers(void)
{
  s=create_int_pointer(l3);
  nv_cw=create_int_pointer(l3);
  nv_co=create_int_pointer(l3);

  nv_w=create_int_pointer(l3);
  nv_o=create_int_pointer(l3);
  nv_s=create_int_pointer(l3);
  nv_g=create_int_pointer(l3);

  inter_pos=create_int_pointer(l3);
  w_inter=create_int_pointer(l3);
  inter_mtx=create_int_pointer(l3);

  return;
}


/*********************************************************************************************
*                                Rotina que libera os ponteiros                              *
*********************************************************************************************/
void free_pointers(void)
{
	free(s);
	free(nv_cw);
	free(nv_co);
	free(nv_w);
	free(nv_o);
	free(nv_s);
	free(nv_g);
	free(inter_pos);
  	free(w_inter);
	free(inter_mtx);

  return;
}

/*********************************************************************************************
*                             Rotina que abre os arquivos de saída                            *
*********************************************************************************************/
void openfiles(void)
{
	
	file_in = 1;

	// Nome do arquivo - o identifier tem que estar correto!
	sprintf(input_file,"%sgota_3d_L_%d_R_%d_Omg_%.1f_w_%d_h_%d_a_%d_LAST.dsf",CI,l,rg,Omg, w,h,a);

 	 // Abrindo entrada
	if((fp_in = fopen(input_file,"r")) == NULL)
	{
      file_in = 0;
    } // Fim do IF

//--------------------------------------------------------------------------------------------
//   IMPLEMENTAR CONTINUIDADE
//--------------------------------------------------------------------------------------------

	if(file_in==1) 
	{
	
		if((fscanf(fp_in,"%c %ld %d",&lixo,&tempo_in, &n_w)!=3)) 
		{
			fprintf(stderr,"\n\n ERRO lendo arq de entrada: primeira linha\n\n");
			exit(-1); 
    	}// Fim do IF


		fprintf(stderr,"Último tempo: %ld\n",tempo_in);
    	numsteps = tempo_in;
    	if(numsteps==mc_steps) 
		{
			fprintf(stderr,"\n\n ESTA AMOSTRA JA ACABOU! SE QUISER AUMENTAR O TEMPO DE \n");
			fprintf(stderr,"       SIMULACAO, AUMENTE A VARIAVEL MC_STEPS!!! \n\n");
			exit(1);
		}// Fim do IF

	    numsteps++;
	    fclose(fp_in);
	    fflush(stdout); 
  
    	// Nome do arquivo de saída 
		sprintf(output_file1,"%sgota_3d_L_%d_R_%d_Omg_%.1f_w_%d_h_%d_a_%d.dsf",CI,l,rg,Omg, w,h,a);
		fp1 = fopen(output_file1,"a");
		fprintf(fp1,"### Continuando..\n");
	    fflush(fp1);
	    fflush(stdout);
	    
		// Arquivo xyz init
        sprintf(output_file1,"%sgota_3d_L_%d_R_%d_Omg_%.1f_w_%d_h_%d_a_%d_init.xyz",CI,l,rg, Omg,w,h,a);	
		fflush(finit);

		// Arquivo xyz
        sprintf(output_file1,"%sgota_3d_L_%d_R_%d_Omg_%.1f_w_%d_h_%d_a_%d.xyz",CI,l,rg,Omg, w,h,a);	
		fxyz = fopen(output_file1,"w");	
		fflush(fxyz); 	   
		
		// Arquivo de Conficguração
		sprintf(output_file1,"%sgota_3d_L_%d_R_%d_Omg_%.1f_w_%d_h_%d_a_%d_conf.dsf", CI,l,rg,Omg, w,h,a);
		fprintf(fconf,"### Continuando..\n");
		fflush(fconf);
		fflush(stdout); 
	    
	    
   
	}// Fim do ELSE

//--------------------------------------------------------------------------------------------
//   Arquivos novos - Para NOVOS runs
//--------------------------------------------------------------------------------------------

	if(file_in==0)  // Se é uma simulação nova
	{
	
		// Arquivo xyz init
        sprintf(output_file1,"%sgota_3d_L_%d_R_%d_Omg_%.1f_w_%d_h_%d_a_%d_init.xyz",CI,l,rg, Omg,w,h,a);	
		finit = fopen(output_file1,"w");	
		fflush(finit);

		// Arquivo xyz
        sprintf(output_file1,"%sgota_3d_L_%d_R_%d_Omg_%.1f_w_%d_h_%d_a_%d.xyz",CI,l,rg,Omg, w,h,a);	
		fxyz = fopen(output_file1,"w");	
		fflush(fxyz);	

      // Arquivo de Conficguração
		sprintf(output_file1,"%sgota_3d_L_%d_R_%d_Omg_%.1f_w_%d_h_%d_a_%d_conf.dsf", CI,l,rg,Omg, w,h,a);
		fconf = fopen(output_file1,"w");

		fprintf(fconf,"# =====================================================================\n");
		fprintf(fconf,"# =                     Parâmetros da Simulação                       =\n");
		fprintf(fconf,"# =====================================================================\n");
		fprintf(fconf,"# Gota - CB (1) ou Wenzel(2)  3D : %d\n",initialstate);
		fprintf(fconf,"#\n");
		fprintf(fconf,"# L = %d   Rg = %d   fo = %3.2f\n" ,l,rg, f_o);
		fprintf(fconf,"# w = %d      h = %d     a = %d\n",w,h,a);
		fprintf(fconf,"#\n");
		fprintf(fconf,"# Omg = %.1f\n" ,Omg);
		fprintf(fconf,"#\n");
		fprintf(fconf,"# Gw = %5.4f   Lambda_w = %5.4f\n",Gw, Lambda_w);
		fprintf(fconf,"# Go = %5.4f   Lambda_o = %5.4f\n",Go, Lambda_o);
		fprintf(fconf,"#\n");
		fprintf(fconf,"# eps_SG = %3.2f     eps_WO = %3.2f\n",eps_SG, eps_WO);
		fprintf(fconf,"# eps_SW = %3.2f     eps_WG = %3.2f\n",eps_SW, eps_WG);
		fprintf(fconf,"# eps_SO = %3.2f     eps_OG = %3.2f\n",eps_SO, eps_OG);
		fprintf(fconf,"#\n");
		fprintf(fconf,"# Temp = %4.2f           kB = %3.2f\n",temp,kB);
		fprintf(fconf,"#\n");
		fprintf(fconf,"# Total time = %d\n",mc_steps);
		fprintf(fconf,"# =====================================================================\n");
		fprintf(fconf,"# t and interface's sites\n");
		fflush(fconf);
		fflush(stdout);		


		// Arquivo dsf - Observáveis  
		sprintf(output_file1,"%sgota_3d_L_%d_R_%d_Omg_%.1f_w_%d_h_%d_a_%d.dsf",CI,l,rg,Omg, w,h,a);
		fflush(stdout);
		fp1 = fopen(output_file1,"w"); //fp1 é o localizador desse arquivo

		fprintf(fp1,"# =====================================================================\n");
		fprintf(fp1,"# =                     Parâmetros da Simulação                       =\n");
		fprintf(fp1,"# =====================================================================\n");
		fprintf(fp1,"# Gota - CB (1) ou Wenzel(2)  3D : %d\n",initialstate);
		fprintf(fp1,"# \n");
		fprintf(fp1,"# L = %d   Rg = %d   fo = %3.2f\n" ,l,rg, f_o);
		fprintf(fp1,"# w = %d      h = %d     a = %d\n",w,h,a);
		fprintf(fp1,"# \n");
		fprintf(fp1,"# Omg = %.1f\n" ,Omg);
		fprintf(fp1,"# \n");
		fprintf(fp1,"# Gw = %5.4f   Lambda_w = %5.4f\n",Gw, Lambda_w);
		fprintf(fp1,"# Go = %5.4f   Lambda_o = %5.4f\n",Go, Lambda_o);
		fprintf(fp1,"# \n");
		fprintf(fp1,"# eps_SG = %3.2f     eps_WO = %3.2f\n",eps_SG, eps_WO);
		fprintf(fp1,"# eps_SW = %3.2f     eps_WG = %3.2f\n",eps_SW, eps_WG);
		fprintf(fp1,"# eps_SO = %3.2f     eps_OG = %3.2f\n",eps_SO, eps_OG);
		fprintf(fp1,"# \n");
		fprintf(fp1,"# Temp = %4.2f           kB = %3.2f\n",temp,kB);
		fprintf(fp1,"# \n");
		fprintf(fp1,"# Total time = %d\n",mc_steps);
		fprintf(fp1,"# =====================================================================\n");
		fprintf(fp1,"# t       V           E    c_x    c_y     Rax   Rrx     Ray   Rry      Bax     Brx    Bay      Bry   Tax   Trx   Tay   Try\n");
   		fflush(fp1);
		fflush(stdout);

  } // Fim do IF

  return;
}
/*********************************************************************************************
*                              Rotina que escreve o arquivo xyz                              *
*********************************************************************************************/
void xyz(int *s, int l, int m)  
{

	int i;
	int l2,l3; 
	int site, npart;
	int x,y,z;

	l2 = l*l;
	l3 = l2*l;

	npart = n_w+n_o+n_s;

	if ( (m == 0) ) 
	{

		fprintf(finit,"%d\n",npart);
		fflush(finit);

		fprintf(finit,"teste\n");
		fflush(finit);

		for(i=0;i<l3;++i)
	    {
			//vizinhos
			site = i;
			x = site%l;
			y = (site/l)%l;
			z = site/l2;

			if ( (s[i] == 1) ) //Água
			{
				fprintf(finit,"O  %d  %d  %d\n",x,y,z);
				fflush(finit);
			}//fim do IF
	
			else if ( (s[i] == 2) ) //óleo
			{
				fprintf(finit,"V  %d  %d  %d\n",x,y,z);
				fflush(finit);
			}//fim do ELSE IF

			else if ( (s[i] == 9 ) ) //óleo
			{
				fprintf(finit,"C  %d  %d  %d\n",x,y,z);
				fflush(finit);
			}//fim do ELSE IF

		}//Fim do for

	}//fim do IF

	else if ( (m == 1) )
	{

		npart = npart;
		fprintf(fxyz,"%d\n",npart);
		fflush(fxyz);

		fprintf(fxyz,"teste\n");
		fflush(fxyz);

		for(i=0;i<l3;++i)
	    {
			//vizinhos
			site = i;
			x = site%l;
			y = (site/l)%l;
			z = site/l2;

			if ( (s[i] == 1) ) //Água
			{
				fprintf(fxyz,"O  %d  %d  %d\n",x,y,z);
				fflush(fxyz);
			}//fim do IF
	
			else if ( (s[i] == 2) ) //óleo
			{
				fprintf(fxyz,"V  %d  %d  %d\n",x,y,z);
				fflush(fxyz);
			}//fim do ELSE IF

			else if ( (s[i] == 9) ) //óleo
			{
				fprintf(fxyz,"C  %d  %d  %d\n",x,y,z);
				fflush(fxyz);
			}//fim do ELSE IF

		}//Fim do for

	}//fim do ELSE IF

 return;

}
/*********************************************************************************************
*                       Rotina que escreve a conf no formato original .dsf                   *
*********************************************************************************************/
void save_conf(int num_steps,int iout) 
{
  
	int i;
  
	if(iout==1) 
	{	
		sprintf(output_file1,"%sgota_3d_L_%d_R_%d_Omg_%.1f_w_%d_h_%d_a_%d_LAST.dsf",CI,l,rg, Omg,w,h,a);
		flast = fopen(output_file1,"w");
		fprintf(flast,"# %d\n",num_steps);


		for(i=0;i<l3;i++) 
		{
			if(s[i]==1 || s[i]==2 ) 
			{
				fprintf(flast,"%d  %d\n",i, s[i]);
			}//fim do IF

		} //fim do FOR

  		fflush(flast);
    
  		fclose(flast);
	} //fim do if

	else if(iout==0) 
	{	

		fprintf(fconf,"# tempo  %d\n",num_steps);
		for(i=0;i<int_label;i++) 
		{
			if(s[w_inter[i]]==1 || s[w_inter[i]]==2 ) 
			{
      			fprintf(fconf,"%d  %d\n",w_inter[i], s[w_inter[i]] );
  			} //fim o IF

		} //fim do FOR

  		fflush(fconf);
  		fprintf(fconf,"\n\n"); 

	} //fim do ELSE IF

	return;
}
/*********************************************************************************************
*                              Rotina que inicializa o programa                              *
*********************************************************************************************/
void initialization(void)
{
	int i;
  
	initial_state(s,l,rg,initialstate);
  
	n_w=0;
	n_o=0;
	n_s=0;
	int_label=0;
	
	for (i=0; i < l3; ++i) 
    {
		if (s[i]==1) ++n_w;
		if (s[i]==2) ++n_o;
		if (s[i]==9) ++n_s;

		if (inter_pos[i]==1) 
		{
			w_inter[int_label]=i;
			inter_mtx[i]=int_label;
			++int_label;

		} // Fim do If

		else inter_mtx[i]=-1;

    } //Fim do FOR

	vol_w = n_w;
	vol_o = n_o;
	vol = vol_w + vol_o;

	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"=                        Início da simulação                        =\n");
	fprintf(fout,"=====================================================================\n");   
	fprintf(fout,"Número de partículas de água    : %d\n",vol_w);
	fprintf(fout,"Número de partículas de óleo    : %d\n",vol_o);
	fprintf(fout,"Número de partículas de líquido : %d\n",vol);
	fprintf(fout,"Número de sítios na interface   : %d\n",int_label);
	fprintf(fout,"=====================================================================\n"); 
	fprintf(fout,"\n"); 
	fflush(stdout);

  return;
}
/*********************************************************************************************
*                          Rotina gera a gota com a distribuição desejada                    *
*********************************************************************************************/
void initial_state(int *s, int L, int Rg, int Initialstate)  
{
  
	int i,j,k;
	int l2,l3,rg2;
	int site;
 	long int count=0;
	int neigh[26];
	int x,y,z,xn,yn,zn;
	double fator;
	double teste;

	l2=L*L;
	l3=l2*L;

	rg2=Rg*Rg; 


	for (i=0;i<l3;++i)
    {
      s[i]=0;
      nv_cw[i]=0; //vizinhos proximos água
      nv_co[i]=0; //vizinhos proximos oleo
      nv_s[i]=0; //vizinhos solido
      nv_w[i]=0; //vizinhos agua
      nv_o[i]=0; //vizinhos oleo
      nv_g[i]=0; //vizinhos gas
      inter_pos[i]=0; 
      w_inter[i]=-1;

    } //Fim do FOR

//--------------------------------------------------------------------------------------------
//   Criando uma nova gota
//--------------------------------------------------------------------------------------------
	if(file_in==0) 
	{

		numsteps=0;
		confstep=0;
		count = 0;

		switch (initialstate)
		{
      	case 1 : // CB

			for (k = 0; k < L; ++k)
				for(j = 0; j < L; ++j)
					for(i = 0; i < L; ++i)
			{

					long z=frac_h*Rg+h+h_base+1;
					site=k*l2+L*j+i;
					if ( (i-L/2+1)*(i-L/2+1) + (j-L/2+1)*(j-L/2+1)+(k-z+1)*(k-z+1) <= rg2)
		  			{

						if ( f_o == 0.0) //Se for apenas água
		  				{
		    				*(s+site)=1;
		    				count++;
						}//fim do IF
						else if ( f_o == 1.0)//Se for apenas óleo
		  				{
		    				*(s+site)=2;
		    				count++;
						}//fim do ELSE IF

						else //se for uma mistura
						{
							/*---------------------------------*/ 
							/*Para gerar uma conf aleatória*/ 

							teste = FRANDOM;
							if ( (teste>f_o)  )
							{

								*(s+site)=1;
		    					count++;

							} //fim do IF
							else if ( (teste<=f_o)  )
							{

		    					*(s+site)=2;
		    					count++;

							} //fim do ELSE IF
							/*---------------------------------*/
							
						}//fim do ELSE

		  			} //fim do IF

			} //Fim do FOR
		break;

		case 2 : // We

			fator = pow(2.0,0.666);
			for (k = 0; k < L; ++k)
				for(j = 0; j < L; ++j)
					for(i = 0; i < L; ++i)
			{

					long z=h_base;
					site=k*l2+l*j+i;
					if ( (i-L/2+1)*(i-L/2+1) + (j-L/2+1)*(j-L/2+1)+(k-z+1)*(k-z+1) <= fator*rg2)
		  			{
		    			if ( f_o == 0.0) //Se for apenas água
		  				{
		    				*(s+site)=1;
		    				count++;
						}//fim do IF
						else if ( f_o == 1.0)//Se for apenas óleo
		  				{
		    				*(s+site)=2;
		    				count++;
						}//fim do ELSE IF

						else //se for uma mistura
						{
							/*---------------------------------*/ 
							/*Para gerar uma conf aleatória*/ 

							teste = FRANDOM;
							if ( (teste>f_o)  )
							{

								*(s+site)=1;
		    					count++;

							} //fim do IF
							else if ( (teste<=f_o)  )
							{

		    					*(s+site)=2;
		    					count++;

							} //fim do ELSE IF
							/*---------------------------------*/
							
						}//fim do ELSE
		  			} //fim do IF

			}//fim do FOR
		break;

		}//Fim do SWITCH

	}//Fim do IF

//--------------------------------------------------------------------------------------------
//   Criando a superfície de pilares
//--------------------------------------------------------------------------------------------

	for (k=0;k<h+h_base;++k) 
	{
		for (j=0;j<L;++j) 
		{
			for (i=0;i<L;++i) 
			{ 
				
				site=k*l2+j*L+i;

				// Superfície de pilares
				if ( (k<h_base) || ((i%(w+a)<w)&&(j%(w+a)<w)) ) s[site]=9;
	
      		} //fim do FOR
		}//fim do FOR   
	} //fim do FOR

//--------------------------------------------------------------------------------------------
//   Criando a lista de vizinhos
//--------------------------------------------------------------------------------------------	   
	for(i=0;i<l3;++i)
    {
		//vizinhos
		site = i;
		x = site%L;
		y = (site/L)%L;
		z = site/l2;
      
		// 0
		yn = (y==0)?L-1:y-1;
		neigh[0]=x + yn*L + z*l2;
		// 1
		yn = (y==L-1)?0:y+1;
		neigh[1] = x + yn*L + z*l2;
		// 2
		xn = (x==L-1)?0:x+1;
		neigh[2] = xn + y*L + z*l2;
		// 3
		xn = (x==0)?L-1:x-1;
		neigh[3] = xn + y*L + z*l2;
		// 4
		zn = (z==0)?L-1:z-1;
		neigh[4] = x + y*L + zn*l2;
		// 5
		zn = (z==L-1)?0:z+1;
		neigh[5] = x + y*L + zn*l2;
		//6 
		xn = (x==L-1)?0:x+1;
		yn = (y==0)?L-1:y-1;
		neigh[6] = xn + yn*L + z*l2;
		// 7
		xn = (x==0)?L-1:x-1;
		yn = (y==0)?L-1:y-1;
		neigh[7] = xn + yn*L + z*l2;
		// 8 
		xn = (x==L-1)?0:x+1;
		yn = (y==L-1)?0:y+1;
		neigh[8] = xn + yn*L + z*l2;
		// 9
		xn = (x==0)?L-1:x-1;
		yn = (y==L-1)?0:y+1;
		neigh[9] = xn + yn*L + z*l2;
		// 10
		xn = (x==L-1)?0:x+1;
		zn = (z==0)?L-1:z-1;
		neigh[10] = xn + y*L + zn*l2;
		// 11
		xn = (x==0)?L-1:x-1;
		zn = (z==0)?L-1:z-1;
		neigh[11] = xn + y*L + zn*l2;
		// 12
		xn = (x==L-1)?0:x+1;
		zn = (z==L-1)?0:z+1;
		neigh[12] = xn + y*L + zn*l2;
		// 13
		xn = (x==0)?L-1:x-1;
		zn = (z==L-1)?0:z+1;
		neigh[13] = xn + y*L + zn*l2;
		// 14
		zn = (z==0)?L-1:z-1;
		yn = (y==0)?L-1:y-1;
		neigh[14] = x + yn*L + zn*l2;
		// 15
		zn = (z==L-1)?0:z+1;
		yn = (y==0)?L-1:y-1;
		neigh[15] = x + yn*L + zn*l2;
		// 16
		zn = (z==0)?L-1:z-1;
		yn = (y==L-1)?0:y+1;
		neigh[16] = x + yn*L + zn*l2;
		// 17
		zn = (z==L-1)?0:z+1;
		yn = (y==L-1)?0:y+1;
		neigh[17] = x + yn*L + zn*l2;
		// 18
		xn = (x==L-1)?0:x+1;
		yn = (y==0)?L-1:y-1;
		zn = (z==0)?L-1:z-1;
		neigh[18] = xn + yn*L + zn*l2;
		// 19
		xn = (x==L-1)?0:x+1;
		yn = (y==0)?L-1:y-1;
		zn = (z==L-1)?0:z+1;
		neigh[19] = xn + yn*L + zn*l2;
		// 20
		xn = (x==0)?L-1:x-1;
		yn = (y==0)?L-1:y-1;
		zn = (z==0)?L-1:z-1;
		neigh[20] = xn + yn*L + zn*l2;
		// 21
		xn = (x==0)?L-1:x-1;
		yn = (y==0)?L-1:y-1;
		zn = (z==L-1)?0:z+1;
		neigh[21] = xn + yn*L + zn*l2;
		// 22
		xn = (x==L-1)?0:x+1;
		yn = (y==L-1)?0:y+1;
		zn = (z==0)?L-1:z-1;
		neigh[22] = xn + yn*L + zn*l2;
		// 23
		xn = (x==L-1)?0:x+1;
		yn = (y==L-1)?0:y+1;
		zn = (z==L-1)?0:z+1;
		neigh[23] = xn + yn*L + zn*l2;
		// 24
		xn = (x==0)?L-1:x-1;
		yn = (y==L-1)?0:y+1;
		zn = (z==0)?L-1:z-1;
		neigh[24] = xn + yn*L + zn*l2;
		// 25	  
		xn = (x==0)?L-1:x-1;
		yn = (y==L-1)?0:y+1;
		zn = (z==L-1)?0:z+1;
		neigh[25] = xn + yn*L + zn*l2;
 
//--------------------------------------------------------------------------------------------
//   Contando o numero de vizinhos de cada tipo
//--------------------------------------------------------------------------------------------	
		for (j=0;j<t_close_neigh;++j)
		{

			if (s[neigh[j]]==1) 
	    	{
				++nv_cw[i];
				++nv_w[i];
	    	} //fim do IF

			else if (s[neigh[j]]==2) 
	    	{
				++nv_co[i];
				++nv_o[i];
	    	} //fim do ELSE IF

			else if (s[neigh[j]]==0) 
	    	{
				++nv_g[i];
	    	}//fim do ELSE IF

			else if (s[neigh[j]]==9)
	    	{
				++nv_s[i];
	    	}//fim do ELSE IF

		}//fim do FOR

//-----------------------------------------------------------------
		for (j=t_close_neigh;j<t_neigh ;++j)
		{
			if (s[neigh[j]]==1) 
	    	{
				++nv_w[i];
	    	} //fim do IF

			else if (s[neigh[j]]==2) 
	    	{
				++nv_o[i];
	    	} //fim do ELSE IF

			else if (s[neigh[j]]==0) 
	    	{
				++nv_g[i];
	    	}//fim do ELSE IF

			else if (s[neigh[j]]==9)
	    	{
				++nv_s[i];
	    	}//fim do ELSE 

		}//fim do FOR

//--------------------------------------------------------------------------------------------
//   Vendo se o sitio esta na interface
//--------------------------------------------------------------------------------------------	
		for (j=0;j<t_close_neigh;++j)
		{

			if ( (s[i]+s[neigh[j]] == 1)  ) //Gas + Água
		    {
				inter_pos[i]=1;
				inter_pos[neigh[j]]=1;
			}//fim do IF
			else if ( (s[i]+s[neigh[j]] == 2) && (s[i]!=s[neigh[j]]) ) //Gas + óleo
		    {	
				inter_pos[i]=1;
				inter_pos[neigh[j]]=1;
			}//fim do ELSE IF
			else if ( (s[i]+s[neigh[j]] == 3) ) //Água + óleo
		    {
				inter_pos[i]=1;
				inter_pos[neigh[j]]=1;
			}//fim do ELSE IF
			else if ( (s[i] == 1)  && (s[neigh[j]] == 9) )  //Água + solido
			{  
				inter_pos[i]=1;
	  		}//fim do ELSE IF
			else if ( (s[i] == 2)  && (s[neigh[j]] == 9)  )  //Oleo + solido
			{  
				inter_pos[i]=1;
	  		}//fim do ELSE IF
	 
		}//fim do FOR
//--------------------------------------------------------------------------------------------

    }//fim do FOR

 return;

}
/*********************************************************************************************
*                               Rotina que faz a dinâmica do MC                              *
*********************************************************************************************/
void dynamics(int *s,int num_steps)
{
	int j;
	int site,s_site, s_teste, label,soma; 
	int xi, hs;
	int nw,ns,no,ng;
	double temp_e, delta_s, delta_v, delta_g, delta_o, delta;
	long int count=0;
	double tg;

	tg = (tan(Omg*M_PI/180.0));

  
	for(j = 0; j < t_vol ; ++j)
    {

//--------------------------------------------------------------------------------------------
//   Sorteando um sítio na interface
//--------------------------------------------------------------------------------------------	

		label = (int) (FRANDOM*int_label);
		site = w_inter[label];
		xi = site%l;
		hs = (site/l2-h_base)+(xi*tg);
		s_site=s[site];

		trial++;
      
		acc_interface++;
      
		acc_teste++;

		delta_s=0;
		delta_v=0;
		delta_o=0;
		delta_g=0;
      
		nw=nv_w[site];
		no=nv_o[site];
		ns=nv_s[site];
		ng=nv_g[site];

//--------------------------------------------------------------------------------------------
//   Determinando qual será a troca
//--------------------------------------------------------------------------------------------	

		if(f_o == 0)
		{

			if(s_site==0)
			{
				s_teste = 1; //troca de ar-agua
			}//fim do IF
			else
			{
				s_teste = 0; //troca de agua-ar
			}//fim do ELSE

		}//fim do IF
		else if(f_o == 1)
		{

			if(s_site==0)
			{
				s_teste = 2; //troca de ar-oleo
			}//fim do IF
			else
			{
				s_teste = 0; //troca de oleo-ar
			}//fim do ELSE

		}//fim do ELSE IF
		else
		{
			s_teste = s_site;
			while (s_teste == s_site)
			{
				s_teste = (int) (FRANDOM*3);
			}//fim do WHILE
			
			if(s_teste==s_site) printf("ERRO: s_teste=s_site!!\n");

		}//fim do ELSE

		soma = s_site +s_teste;	

//--------------------------------------------------------------------------------------------
//   Calculando a energia da troca
//--------------------------------------------------------------------------------------------
		
		if (soma==1)  //água e ar
		{

			if (s_site==0)
			{

				delta_g = Gw*hs;
				delta_s = (ng-nw)*eps_WG + ns*(eps_SW-eps_SG) + (eps_WO-eps_OG)*no;
				delta_v = (1+2*(vol_w-t_vol_w))*Lambda_w; // Ganha uma água

			}//fim do IF

			else //if (s_site==1)
			{

				delta_g = -Gw*hs; 
				delta_s = (nw-ng)*eps_WG + ns*(eps_SG-eps_SW) + (eps_OG-eps_WO)*no;
				delta_v = (1-2*(vol_w-t_vol_w))*Lambda_w; // Perco uma água

			}// fim do ELSE
			
		} //fim do IF

		if ((soma==2) && (s_teste!=s_site))  //óleo e ar
		{

			if (s_site==0)
			{

				delta_g = Go*hs;
				delta_s = (ng-no)*eps_OG + ns*(eps_SO-eps_SG) + (eps_WO-eps_WG)*nw;
				delta_o = (1+2*(vol_o-t_vol_o))*(Lambda_o); // Ganho um um óleo
	  

			}//fim do IF

			else //if (s_site==2)
			{
	
				delta_g = -Go*hs; 
				delta_s = (no-ng)*eps_OG + ns*(eps_SG-eps_SO) + (eps_WG-eps_WO)*nw;
				delta_o = (1-2*(vol_o-t_vol_o))*Lambda_o;  // Perco um  um óleo

						  

			}// fim do ELSE

		} //fim do  IF

		if (soma==3)  //óleo e água
		{
			
			if (s_site==1)
			{	
				delta_g = (Go-Gw)*hs;
				delta_s = (nw-no)*eps_WO + (eps_OG-eps_WG)*ng + (eps_SO-eps_SW)*ns;
				delta_o = (1+2*(vol_o-t_vol_o))*Lambda_o; //Ganho um óleo
				delta_v = (1-2*(vol_w-t_vol_w))*Lambda_w; // Perco uma água
				

			}//fim do IF

			else //if (s_site==2)
			{
				delta_g = (Gw-Go)*hs; 
				delta_s = (no-nw)*eps_WO + (eps_WG-eps_OG)*ng + (eps_SW-eps_SO)*ns;
				delta_o = (1-2*(vol_o-t_vol_o))*Lambda_o;  // Perco um óleo
				delta_v = (1+2*(vol_w-t_vol_w))*Lambda_w; // Ganha uma água

			}// fim do ELSE	

		}//fim do IF

		if (soma>3)
		{
			printf("ENTROU!!\n");
			delta_g = 100000.0; 
			delta_s = 100000.0;
			delta_v = 100000.0;
		}//fim do else
 
		delta=delta_v+delta_s+delta_g+delta_o;

//--------------------------------------------------------------------------------------------
//  Testando a troca
//--------------------------------------------------------------------------------------------
      
		if (delta<=0) 
		{

			flip_spin(site,s_teste);
			count++;
			accepted++;

		}//fim do IF

 		else 
		{
			if (T > 0.0) 
			{
				temp_e=FRANDOM;
				if(temp_e<exp(-delta/(kB*T))) 
				{

					flip_spin(site,s_teste);
					count++;
					accepted++;

				}//fim do IF
			}//fim do IF
		}//dim do ELSE

	} //fim do FOR

  return;
}
/*********************************************************************************************
*                               Rotina que flipa o spin do sítio                             *
*********************************************************************************************/
void flip_spin (long site, long new_s)
{
	int i,neighbour;
	int soma;
	int neigh[26];
	int x,y,z,xn,yn,zn;

//--------------------------------------------------------------------------------------------
//  Atualizando o volume
//--------------------------------------------------------------------------------------------

	soma = s[site]+new_s;

	if (soma==1)  //água e ar
	{
		if (s[site]==0)
		{
			++vol;
			++vol_w;
		}//fim do IF

		else //if (s_site==1)
		{
			--vol;
			--vol_w;
		}// fim do ELSE		
	} //fim do IF

	else if (soma==2)  //óleo e ar
	{
		if (s[site]==0)
		{
			++vol;
			++vol_o;
		}//fim do IF

		else //if (s_site==1)
		{
			--vol;
			--vol_o;
		}// fim do ELSE		

	} //fim do ELSE IF

	else if (soma==3)  //óleo e água
	{
		if (s[site]==1)
		{
			++vol_o;
			--vol_w;
		}//fim do IF

		else //if (s_site==2)
		{
			--vol_o;
			++vol_w;
		}// fim do ELSE					
					
	}//fim do ELSE IF

	n_w = (int) vol_w;
	n_o = (int) vol_o;

	
//--------------------------------------------------------------------------------------------
//  Atualizando o numero de vizinhos
//--------------------------------------------------------------------------------------------
  
	x = site%l;
	y = (site/l)%l;
	z = site/l2;
      
	// 0
	yn = (y==0)?l-1:y-1;
	neigh[0]=x + yn*l + z*l2;
	// 1
	yn = (y==l-1)?0:y+1;
	neigh[1] = x + yn*l + z*l2;
	// 2
	xn = (x==l-1)?0:x+1;
	neigh[2] = xn + y*l + z*l2;
	// 3
	xn = (x==0)?l-1:x-1;
	neigh[3] = xn + y*l + z*l2;
	// 4
	zn = (z==0)?l-1:z-1;
	neigh[4] = x + y*l + zn*l2;
	// 5
	zn = (z==l-1)?0:z+1;
	neigh[5] = x + y*l + zn*l2;
	//6 
	xn = (x==l-1)?0:x+1;
	yn = (y==0)?l-1:y-1;
	neigh[6] = xn + yn*l + z*l2;
	// 7
	xn = (x==0)?l-1:x-1;
	yn = (y==0)?l-1:y-1;
	neigh[7] = xn + yn*l + z*l2;
	// 8 
	xn = (x==l-1)?0:x+1;
	yn = (y==l-1)?0:y+1;
	neigh[8] = xn + yn*l + z*l2;
	// 9
	xn = (x==0)?l-1:x-1;
	yn = (y==l-1)?0:y+1;
	neigh[9] = xn + yn*l + z*l2;
	// 10
	xn = (x==l-1)?0:x+1;
	zn = (z==0)?l-1:z-1;
	neigh[10] = xn + y*l + zn*l2;
	// 11
	xn = (x==0)?l-1:x-1;
	zn = (z==0)?l-1:z-1;
	neigh[11] = xn + y*l + zn*l2;
	// 12
	xn = (x==l-1)?0:x+1;
	zn = (z==l-1)?0:z+1;
	neigh[12] = xn + y*l + zn*l2;
	// 13
	xn = (x==0)?l-1:x-1;
	zn = (z==l-1)?0:z+1;
	neigh[13] = xn + y*l + zn*l2;
	// 14
	zn = (z==0)?l-1:z-1;
	yn = (y==0)?l-1:y-1;
	neigh[14] = x + yn*l + zn*l2;
	// 15
	zn = (z==l-1)?0:z+1;
	yn = (y==0)?l-1:y-1;
	neigh[15] = x + yn*l + zn*l2;
	// 16
	zn = (z==0)?l-1:z-1;
	yn = (y==l-1)?0:y+1;
	neigh[16] = x + yn*l + zn*l2;
	// 17
	zn = (z==l-1)?0:z+1;
	yn = (y==l-1)?0:y+1;
	neigh[17] = x + yn*l + zn*l2;
	// 18
	xn = (x==l-1)?0:x+1;
	yn = (y==0)?l-1:y-1;
	zn = (z==0)?l-1:z-1;
	neigh[18] = xn + yn*l + zn*l2;
	// 19
	xn = (x==l-1)?0:x+1;
	yn = (y==0)?l-1:y-1;
	zn = (z==l-1)?0:z+1;
	neigh[19] = xn + yn*l + zn*l2;
	// 20
	xn = (x==0)?l-1:x-1;
	yn = (y==0)?l-1:y-1;
	zn = (z==0)?l-1:z-1;
	neigh[20] = xn + yn*l + zn*l2;
	// 21
	xn = (x==0)?l-1:x-1;
	yn = (y==0)?l-1:y-1;
	zn = (z==l-1)?0:z+1;
	neigh[21] = xn + yn*l + zn*l2;
	// 22
	xn = (x==l-1)?0:x+1;
	yn = (y==l-1)?0:y+1;
	zn = (z==0)?l-1:z-1;
	neigh[22] = xn + yn*l + zn*l2;
	// 23
	xn = (x==l-1)?0:x+1;
	yn = (y==l-1)?0:y+1;
	zn = (z==l-1)?0:z+1;
	neigh[23] = xn + yn*l + zn*l2;
	// 24
	xn = (x==0)?l-1:x-1;
	yn = (y==l-1)?0:y+1;
	zn = (z==0)?l-1:z-1;
	neigh[24] = xn + yn*l + zn*l2;
	// 25	  
	xn = (x==0)?l-1:x-1;
	yn = (y==l-1)?0:y+1;
	zn = (z==l-1)?0:z+1;
	neigh[25] = xn + yn*l + zn*l2;

	for(i=0;i<t_close_neigh;++i)
    {
		neighbour=neigh[i];

		if (soma==1)
		{
			if (new_s==1)
			{
				++nv_w[neighbour];
				++nv_cw[neighbour];
				--nv_g[neighbour];
			} //fim do IF
			else 
			{
				--nv_w[neighbour];
				--nv_cw[neighbour];
				++nv_g[neighbour];
			}//fim do ELSE
		}//fim do IF

		else if (soma==2)
		{
			if (new_s==2)
			{
				++nv_o[neighbour];
				++nv_co[neighbour];
				--nv_g[neighbour];
			} //fim do IF
			else 
			{
				--nv_o[neighbour];
				--nv_co[neighbour];
				++nv_g[neighbour];
			}//fim do ELSE
		}//fim do ELSE IF

		else if (soma==3)
		{
			if (new_s==1)
			{
				++nv_w[neighbour];
				++nv_cw[neighbour];
				--nv_o[neighbour];
				--nv_co[neighbour];

			} //fim do IF
			else 
			{
				--nv_w[neighbour];
				--nv_cw[neighbour];
				++nv_o[neighbour];
				++nv_co[neighbour];
			}//fim do ELSE
		}//fim do ELSE IF

	}//fim do FOR


	for (i=t_close_neigh;i<t_neigh ;++i)
    {
		neighbour=neigh[i];

		if (soma==1)
		{
			if (new_s==1)
			{
				++nv_w[neighbour];
				--nv_g[neighbour];
			} //fim do IF
			else 
			{
				--nv_w[neighbour];
				++nv_g[neighbour];
			}//fim do ELSE
		}//fim do IF

		else if (soma==2)
		{
			if (new_s==2)
			{
				++nv_o[neighbour];
				--nv_g[neighbour];
			} //fim do IF
			else 
			{
				--nv_o[neighbour];
				++nv_g[neighbour];
			}//fim do ELSE
		}//fim do ELSE IF

		else if (soma==3)
		{
			if (new_s==1)
			{
				++nv_w[neighbour];
				--nv_o[neighbour];

			} //fim do IF
			else 
			{
				--nv_w[neighbour];
				++nv_o[neighbour];
			}//fim do ELSE
		}//fim do ELSE IF

    }//fim do FOR

//--------------------------------------------------------------------------------------------
//  Atualizando sítios na interface
//--------------------------------------------------------------------------------------------

	for (i=0;i<t_close_neigh;++i)
	{
		neighbour=neigh[i];

		if (new_s==0)
		{
			if ( (s[neighbour]==1) || (s[neighbour]==2) )
			{
				if (inter_pos[neighbour]==0) {add_to_interface(neighbour);}
				if (inter_pos[site]==0) {add_to_interface(site);}

	    	}//fim do IF

			else if ((s[neighbour]==0))
			{
	   			if ((nv_cw[neighbour]==0) && (nv_co[neighbour]==0) && (inter_pos[neighbour]==1))
				{
					remove_from_interface(neighbour);
				} //fim do IF
			}//fim do ELSE IF

		}//fim do IF
		
		else if(new_s==1)
		{
			if (s[neighbour]==0 || (s[neighbour]==2))
			{
				if (inter_pos[neighbour]==0) {add_to_interface(neighbour);}
				if (inter_pos[site]==0) {add_to_interface(site);}
			}//fim do IF

			else if (s[neighbour]==1)
			{

			    if ((nv_cw[neighbour]==t_close_neigh) && (inter_pos[neighbour]==1) )
				{
					remove_from_interface(neighbour);
				} //fim do IF
			}//fim do ELSE IF
		}//fim do ELSE IF

		else if(new_s==2)
		{
			if (s[neighbour]==0 || (s[neighbour]==1))
			{
				if (inter_pos[neighbour]==0) {add_to_interface(neighbour);}
				if (inter_pos[site]==0) {add_to_interface(site);}
			}//fim do IF

			else if (s[neighbour]==2)
			{

			    if ((nv_co[neighbour]==t_close_neigh) && (inter_pos[neighbour]==1) )
				{
					remove_from_interface(neighbour);
				} //fim do IF
			}//fim do ELSE IF
		}//fim do ELSE IF

	}//fim do FOR
      
	if (new_s==1) 
	{
		if ((nv_cw[site]==t_close_neigh)&&(inter_pos[site]==1)) {remove_from_interface(site);}
    }//fim do IF
	else if (new_s==2)
	{
		if ((nv_co[site]==t_close_neigh)&&(inter_pos[site]==1)) {remove_from_interface(site);}
	}//fim do ELSE IF
  	else if (new_s==0)
    {
      if (((nv_co[site]+nv_cw[site])==0)&&(inter_pos[site]==1)) {remove_from_interface(site);}
    }//fim do ELSE IF
  
	s[site]=new_s;
     
  return;
}
/*********************************************************************************************
*                      Rotina que adiciona um sítio da lista da interface                    *
*********************************************************************************************/
void add_to_interface(int site)
{
 
	inter_pos[site]=1;
	inter_mtx[site]=int_label;
	w_inter[int_label]=site;
	++int_label;

  return;
}
/*********************************************************************************************
*                       Rotina que remove um sítio da lista da interface                     *
*********************************************************************************************/
void remove_from_interface(int site)
{
	int site_label, last_interface_site;
  
	inter_pos[site]=0;

  	site_label=inter_mtx[site];
  	inter_mtx[site]=-1;
  
	last_interface_site=w_inter[int_label-1];
  	w_inter[site_label]=last_interface_site;
	inter_mtx[last_interface_site]=site_label;
  
	w_inter[int_label-1]=-1;
	--int_label;
  
  return;

}
/*********************************************************************************************
*                                Rotina que calcula observáveis                              *
*********************************************************************************************/
void measure_angle(int num_steps) 
{

	int i,j;
	int x,y,z;
	int site;
	int xmin,xmax,ymin,ymax;	

// -------------------------------------------------------------------------------------------
// Determinando a altura da gota


	drop_height=0;
	drop_bottom=l;
	center_x=0;
	center_y=0;
	
	for(i=0;i<int_label;++i)
    {

		site = w_inter[i];
		if ( s[site]==1 )
		{
			x = site%l;
			y = (site/l)%l;
			z = (site/l2);

			if (z>drop_height) 
			{

				drop_height=z;
				center_x=x;
				center_y=y;

			}//fim do IF
			if (z<drop_bottom) drop_bottom=z;

		}//fim do IF

	}//Fim do FOR

// -------------------------------------------------------------------------------------------
// Determinando o xmax, xmin, ymax, ymin

	xmin=l;
	xmax=0;
	ymin=l;
	ymax=0;

	for(i=0;i<l;++i)
		for(j=0;j<l;++j)
	{

		site = (h_base+h)*l2+j*l+i;

		if (s[site]==1)
		{
			if (i<xmin) xmin=i;
			if (i>xmax) xmax=i;
			if (j<ymin) ymin=j;
			if (j>ymax) ymax=j;

		} //fim do IF

	} //fim dos FOR

// -------------------------------------------------------------------------------------------
// Determinando os raios da base

	Bax = center_x - xmin;
	if(Bax<0.0) Bax=0.0;
	if(Bax>l/2) Bax=(float)l/2;

	Bay = center_y - ymin;
	if(Bay<0.0) Bay=0.0;
	if(Bay>l/2) Bay=(float)l/2;

	Brx = xmax - center_x;
	if(Brx<0.0) Brx=0.0;
	if(Brx>l/2) Brx=(float)l/2;

	Bry = ymax - center_y;
	if(Bry<0.0) Bry=0.0;
	if(Bry>l/2) Bry=(float)l/2;

// -------------------------------------------------------------------------------------------
// Calculando as hipotenusas

	Hax = drop_height*drop_height + Bax*Bax;
	Hay = drop_height*drop_height + Bay*Bay;

	Hrx = drop_height*drop_height + Brx*Brx;
	Hry = drop_height*drop_height + Bry*Bry;
	
// -------------------------------------------------------------------------------------------
// Calculando os raios

	Rax = Hax/(2.0*drop_height);
	if(Rax>(l/2)) Rax=(float)l/2;

	Ray = Hay/(2.0*drop_height);
	if(Ray>(l/2)) Ray=(float)l/2;

	Rrx = Hrx/(2.0*drop_height);
	if(Rrx>(l/2)) Rrx=(float)l/2;

	Rry = Hry/(2.0*drop_height);
	if(Rry>(l/2)) Rry=(float)l/2;

// -------------------------------------------------------------------------------------------
// Calculando os ângulos

	
	Tax = asin(Bax/Rax);
	Tay = asin(Bay/Ray);

	Trx = asin(Brx/Rrx);
	Try = asin(Bry/Rry);

	if (drop_height>Rax) Tax = M_PI - Tax;
	if (drop_height>Ray) Tay = M_PI - Tay;

	if (drop_height>Rrx) Trx = M_PI - Trx;
	if (drop_height>Rry) Try = M_PI - Try;

	Tax = Tax*180.0/M_PI; 
	Tay = Tay*180.0/M_PI;

	Trx = Trx*180.0/M_PI; 
	Try = Try*180.0/M_PI;


// -------------------------------------------------------------------------------------------
  
	fprintf(fp1,"%8d %d %f %d %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n", num_steps, (int)vol, calculate_energy(num_steps),center_x,center_y, Rax,Rrx,Ray,Rry, Bax,Brx,Bay,Bry, Tax,Trx,Tay,Try);
	fflush(fp1);

  return;

}

/*********************************************************************************************
*                                 Rotina que calcula a energia                               *
*********************************************************************************************/
double calculate_energy(int num_steps)
{
	double total=0.0;
	int i,j;
	int hs;
	int neigh[26];  
	int x,y,z,xn,yn,zn,site;
	long int num_interface=0;
	double ene_interface=0.0,ene_grav=0.0;
	double tg;

	tg = (tan(Omg*M_PI/180.0));
	
	for(i=0;i<l3;++i)
    {
		//vizinhos
		site = i;
		x = site%l;
		y = (site/l)%l;
		z = site/l2;
		hs = x*tg;
      
		// 0
		yn = (y==0)?l-1:y-1;
		neigh[0]=x + yn*l + z*l2;
		// 1
		yn = (y==l-1)?0:y+1;
		neigh[1] = x + yn*l + z*l2;
		// 2
		xn = (x==l-1)?0:x+1;
		neigh[2] = xn + y*l + z*l2;
		// 3
		xn = (x==0)?l-1:x-1;
		neigh[3] = xn + y*l + z*l2;
		// 4
		zn = (z==0)?l-1:z-1;
		neigh[4] = x + y*l + zn*l2;
		// 5
		zn = (z==l-1)?0:z+1;
		neigh[5] = x + y*l + zn*l2;
		//6 
		xn = (x==l-1)?0:x+1;
		yn = (y==0)?l-1:y-1;
		neigh[6] = xn + yn*l + z*l2;
		// 7
		xn = (x==0)?l-1:x-1;
		yn = (y==0)?l-1:y-1;
		neigh[7] = xn + yn*l + z*l2;
		// 8 
		xn = (x==l-1)?0:x+1;
		yn = (y==l-1)?0:y+1;
		neigh[8] = xn + yn*l + z*l2;
		// 9
		xn = (x==0)?l-1:x-1;
		yn = (y==l-1)?0:y+1;
		neigh[9] = xn + yn*l + z*l2;
		// 10
		xn = (x==l-1)?0:x+1;
		zn = (z==0)?l-1:z-1;
		neigh[10] = xn + y*l + zn*l2;
		// 11
		xn = (x==0)?l-1:x-1;
		zn = (z==0)?l-1:z-1;
		neigh[11] = xn + y*l + zn*l2;
		// 12
		xn = (x==l-1)?0:x+1;
		zn = (z==l-1)?0:z+1;
		neigh[12] = xn + y*l + zn*l2;
		// 13
		xn = (x==0)?l-1:x-1;
		zn = (z==l-1)?0:z+1;
		neigh[13] = xn + y*l + zn*l2;
		// 14
		zn = (z==0)?l-1:z-1;
		yn = (y==0)?l-1:y-1;
		neigh[14] = x + yn*l + zn*l2;
		// 15
		zn = (z==l-1)?0:z+1;
		yn = (y==0)?l-1:y-1;
		neigh[15] = x + yn*l + zn*l2;
		// 16
		zn = (z==0)?l-1:z-1;
		yn = (y==l-1)?0:y+1;
		neigh[16] = x + yn*l + zn*l2;
		// 17
		zn = (z==l-1)?0:z+1;
		yn = (y==l-1)?0:y+1;
		neigh[17] = x + yn*l + zn*l2;
		// 18
		xn = (x==l-1)?0:x+1;
		yn = (y==0)?l-1:y-1;
		zn = (z==0)?l-1:z-1;
		neigh[18] = xn + yn*l + zn*l2;
		// 19
		xn = (x==l-1)?0:x+1;
		yn = (y==0)?l-1:y-1;
		zn = (z==l-1)?0:z+1;
		neigh[19] = xn + yn*l + zn*l2;
		// 20
		xn = (x==0)?l-1:x-1;
		yn = (y==0)?l-1:y-1;
		zn = (z==0)?l-1:z-1;
		neigh[20] = xn + yn*l + zn*l2;
		// 21
		xn = (x==0)?l-1:x-1;
		yn = (y==0)?l-1:y-1;
		zn = (z==l-1)?0:z+1;
		neigh[21] = xn + yn*l + zn*l2;
		// 22
		xn = (x==l-1)?0:x+1;
		yn = (y==l-1)?0:y+1;
		zn = (z==0)?l-1:z-1;
		neigh[22] = xn + yn*l + zn*l2;
		// 23
		xn = (x==l-1)?0:x+1;
		yn = (y==l-1)?0:y+1;
		zn = (z==l-1)?0:z+1;
		neigh[23] = xn + yn*l + zn*l2;
		// 24
		xn = (x==0)?l-1:x-1;
		yn = (y==l-1)?0:y+1;
		zn = (z==0)?l-1:z-1;
		neigh[24] = xn + yn*l + zn*l2;
		// 25	  
		xn = (x==0)?l-1:x-1;
		yn = (y==l-1)?0:y+1;
		zn = (z==l-1)?0:z+1;
		neigh[25] = xn + yn*l + zn*l2;

    
		for(j=0;j<t_neigh;++j)
		{
			if (s[i]!=s[neigh[j]])
			{
				num_interface++;

				switch (s[i]+s[neigh[j]])
				{
				// s=1 (agua) -- s=0 (ar) -- s=2 (oleo)  --- s=9 (solido)
				case 1: //agua - ar
				total+=0.5*(eps_WG); 
				ene_interface+=0.5*(eps_WG);
				break;
				case 2: //oleo -ar
				total+=0.5*(eps_OG);
				ene_interface+=0.5*(eps_OG);
				break;
				case 3: //agua - oleo
				total+=0.5*(eps_WO);
				ene_interface+=0.5*(eps_WO);
				break;
				case 9: //ar - solido
				total+=0.5*(eps_SG);
				ene_interface+=0.5*(eps_SG);
				break;
				case 10: //agua - solido
				total+=0.5*(eps_SW);
				ene_interface+=0.5*(eps_SW);
				break;
				case 11: //solido - oleo
				total+=0.5*(eps_SO);
				ene_interface+=0.5*(eps_SO);
				break;
				} //fim do SWITCH

			}//fim do IF

		}//fim do FOR


		if (s[i]==1) 
		{
			total+=Gw*(z+hs);
			ene_grav+=Gw*(z+hs);
		}//fim do if

		if (s[i]==2) 
		{
			total+=Go*(z+hs);
			ene_grav+=Go*(z+hs);
		}//fim do if


	} //fim do FOR

	total+=Lambda_w*(vol-t_vol)*(vol-t_vol);
	total+=Lambda_o*(vol_o-t_vol_o)*(vol_o-t_vol_o);


	fprintf(fout,"Total          : %f \n",total);
 	fprintf(fout,"Interfacial    : %f        num_interface: %ld\n",ene_interface,num_interface/2/t_neigh );
	fprintf(fout,"Gravitacional  : %f             int/grav: %f\n",ene_grav,ene_interface/ene_grav);
	fprintf(fout,"Lagrange Total : %f \n",Lambda_w*(vol-t_vol)*(vol-t_vol));
	fprintf(fout,"Lagrange Fração: %f \n",Lambda_o*(vol_o-t_vol_o)*(vol_o-t_vol_o));
	fprintf(fout,"\n\n\n");

	return total;
}

