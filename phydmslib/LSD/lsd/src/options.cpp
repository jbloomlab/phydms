//    LSD - Least Square Dating for etimating substitution rate and divergence dates
//    Copyright (C) <2015> <Thu-Hien To>

//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#include "options.h"

options* getOptions( int argc, char** argv  , bool & rooted)
{
	if( argc>1 )
		return getCommandLine( argc, argv , rooted);
	else
		return getInterface(rooted);
}

options* getCommandLine( int argc, char** argv , bool& rooted)
{
	options* opt = (struct options *) malloc( sizeof(struct options) );
	setDefaultOptions( opt );
	int c;
    string s,path1;
	bool iflag = false, 
	     dflag = false,
         flagA=false,
         flagZ=false;
    char* fnOut;
    while ( (c = getopt(argc, argv, ":i:d:o:s:n:g:r:vct:w:b:ha:z:")) != -1 ) 
    {
        switch (c) 
        {
			case 'i':
				if( access( optarg, R_OK )!=0 )
					myExit( "Cannot read the file named \"%s\"\n", optarg );
				opt->inFile = optarg;	
				iflag = true;
				break;
			case 'd':
                      if( access( optarg, R_OK )!=0 )
                          myExit( "Cannot read the file named \"%s\"\n", optarg );
                      opt->inDateFile = optarg;	
                 dflag = true;
				break;
			case 'o':
                fnOut =new char[strlen(optarg)];
                strcpy(fnOut,optarg);
                opt->outFile = new char[strlen(fnOut)];
                strcpy(opt->outFile,fnOut);
                delete[] fnOut;
				break;
			case 's':
				if( !isInteger(optarg) )
					myExit("Argument of option -s must be an integer.\n");
				opt->seqLength = atoi(optarg);
				if( opt->seqLength<1 )
					myExit("Argument of option -s must be strictly positive.\n");
				break;
			case 'n':
				if( !isInteger(optarg) )
					myExit("Argument of option -n must be an integer.\n");
				opt->nbData = atoi(optarg);
				if( opt->nbData<1 )
					myExit("Argument of option -n must be strictly positive.\n");
				break;
			case 'g':
			     if( access( optarg, R_OK )!=0 )
				     myExit( "Cannot read the file named \"%s\"\n", optarg );
                 opt->fnOutgroup = optarg;  
                 //opt->inFileOp = getDefaultInFileOpName(opt->inFile);
				 break;
		    case 'r':
                 opt->estimate_root = optarg;
                 if (strcmp(opt->estimate_root,"l")!=0 && strcmp(opt->estimate_root,"a")!=0 && strcmp(opt->estimate_root,"as")!=0)
                    myExit("Argument of option -r must be either \"l\" or \"a\" or \"as\".\n");
                 break;
			case 'v':
				opt->variance = true;
				break;
			case 'c':
                 opt->constraint = true;
                 break;
            case 'b':
                 opt->c = atoi( optarg );     
                 if (opt->delta<0)
                    myExit("Argument of option -b must not be negative.\n");
                 break;
            case 't':
                 opt->delta = atof(optarg);
                 if (opt->delta<=0)
                    myExit("Argument of option -t must be strictly positive.\n");
                 break;
            case 'w':
			     if( access( optarg, R_OK )!=0 )
				     myExit( "Cannot read the file named \"%s\"\n", optarg );
                 opt->rate = optarg; 
				 break;
			case 'a':
                 opt->mrca=atof(optarg);
                 flagA=true;
                 break;
            case 'z':
                 opt->leaves=atof(optarg);
                 flagZ=true;
                 break;
			case 'h':
				printHelp();
				exit( EXIT_SUCCESS );				
			case '?':
				myExit("Unrecognized option: -%c\n", optopt);
            case ':':
				myExit("Option -%c requires an operand\n", optopt );
			default:
				myExit("?? getopt returned character code 0%o ??\n", c );
        }
    }
    if( !(iflag) )
		myExit("Argument -i is necessary to continue...\n");
	if (!dflag && (flagA && !flagZ)){
       myExit("The input date file is not provided, so option -z is needed to use with option -a to estimate relative dates.\n");
    }		
	if (!dflag && (!flagA && flagZ)){
       myExit("The input date file is not provided, so option -a is needed to use with option -z to estimate relative dates.\n");
    }		    
	if (!dflag && (!flagA && !flagZ)){
       opt->relative=true;
       opt->mrca=0;
       opt->leaves=1;         
    }
    if (!dflag && (flagA && flagZ)){
       if (opt->mrca >= opt->leaves) 
          myExit("The root date must be strictly smaller than the tips date.\n");
	   opt->relative=true;
    }
    if (dflag) 
       opt->relative=false;   
    int n,m;
    counting(opt->inFile,n,m,1);
    if (m!=2*n && m!=2*n+1) {
        cout<<"Errors in the input trees"<<endl;
        exit(EXIT_FAILURE);
    }
    if (m==2*n+1){
        if (opt->estimate_root==NULL && opt->fnOutgroup==NULL) {
            cout<<"The input trees are not rooted, use either option -g to specify the outgroups file or -r to estimate the root"<<endl;
            exit(EXIT_FAILURE);
        }
        else{
            rooted=false;
        }
    }
    if (opt->fnOutgroup!=NULL) {
        rooted=true;
    }
    if (!rooted && opt->estimate_root!=NULL && strcmp(opt->estimate_root,"l")==0) {
       cout<<"The input trees are not rooted, so the \"l\" method for rooting function is turned to the \"a\" method."<<endl;
       strcpy(opt->estimate_root,"a");
    }
    if (!opt->constraint && opt->estimate_root!=NULL && strcmp(opt->estimate_root,"as")==0){
       cout<<"The non constrainted mode is chosen, so the \"as\" method for rooting function is the same as the \"a\" method."<<endl;                    
       strcpy(opt->estimate_root,"a");  
    }    
    if( opt->outFile==NULL ) opt->outFile = getDefaultOutputFileName( (string)(opt->inFile) );
    return opt;	
}

char* getDefaultOutputFileName( string in)
{
    char *out = (char *) malloc( sizeof(char)*(in.length()+12) );
	strcpy( out, in.c_str() );
	strcat( out, ".result" );
	return out;
}

char* getDefaultOutputNexusTreeFileName( string in)
{
    char *out = (char *) malloc( sizeof(char)*(in.length()+11) );
	strcpy( out, in.c_str()	);
	strcat( out, ".nexus" );
	return out;
}

char* getDefaultOutputNewick1TreeFileName( string in)
{
    char *out = (char *) malloc( sizeof(char)*(in.length()+16) );
	strcpy( out, in.c_str() );
	strcat( out, ".newick" );
	return out;
}

char* getDefaultOutputNewick2TreeFileName( string in)
{
    char *out = (char *) malloc( sizeof(char)*(in.length()+17) );
	strcpy( out, in.c_str() );
	strcat( out, ".date.newick" );
	return out;
}

char * getDefaultInFileOpName( string in){
     char *out = (char *) malloc( sizeof(char)*(in.length()+13) );;
     strcpy(out,in.c_str());
     strcat(out,".ingroup");
     return out;     
}

options* getInterface(bool & rooted)
{
	options* opt = (struct options *) malloc( sizeof(struct options) );
    setDefaultOptions( opt );
    bool inputrooted=true;
	opt->inFile = getInputTreeFileName("Enter your Input Tree File name> ",inputrooted);
    if (!inputrooted){
       opt->estimate_root=new char[2];
       strcpy(opt->estimate_root,"a");
       rooted=false;
    }
    opt->inDateFile = getInputFileName("Enter your Input Date File name if there is any, otherwise press Enter> ");
    cout<<endl;
    if (string(opt->inDateFile).length()==0){               
       cout<<"There is no date file, so the program will estimate relative dates with root date = 0 and tips date = 1."<<endl;
	   cout<<"Type 'y' to continue or 'n' to modify the root date and the tips date"<<endl;		
	   char letter[3];
	   do {
		  fgets( letter, 3, stdin );	   
		  if (*letter=='n' || *letter=='N'){
			do {           
				opt->mrca = getInputReal("Enter the root date (default=0)> ");
				opt->leaves = getInputReal("Enter the tips date (default=1)> ");        
				if (opt->leaves <= opt->mrca) cout<<"Root date must be smaller than the tips date."<<endl;
			} while (opt->leaves <= opt->mrca);
		  }
		  else if (*letter=='y' || *letter=='Y'){
				opt->mrca=0;
				opt->leaves=1;
		 }
		 else {
			cout<<"Type 'y' to continue or 'n' to modify the root date and tips date"<<endl;		
		}		
	   } while (*letter!='n' && *letter!='N' && *letter!='y' && *letter!='Y');
	   opt->relative=true;
    }
    else opt->relative=false;
    opt->outFile = getDefaultOutputFileName( (string)(opt->inFile) );
    bool out=false;
	char letter[3];
	do
	{
		printInterface( stdout, opt , rooted);
		cout<<endl;
		fgets( letter, 3, stdin );
		if( isOptionActivate( opt, *letter ) )
			setOptionsWithLetter( opt, *letter , inputrooted , out);
        rooted = inputrooted || out;
	} while( *letter!='y' && *letter!='Y' );
	return opt;
}


void printInterface( FILE* in, options* opt , bool rooted)
{
    fprintf(in,"\nLEAST-SQUARE METHODS TO ESTIMATE RATES AND DATES - " VERSION" \n\n");
	fprintf(in,"\nInput files:\n");
	fprintf(in,"  I                                 Input tree file : %s\n",opt->inFile);
	if (opt->relative==true)
       fprintf(in,"  D                         Estimate relative dates : mrca date = %.3f, tips date =%.3f\n",opt->mrca,opt->leaves);
	else
	    fprintf(in,"  D                                 Input date file : %s\n",opt->inDateFile);
	fprintf(in,"Output file:\n");	 
	fprintf(in,"  O                                    Output file  : %s\n",opt->outFile);
	fprintf(in,"Parameters:\n");
	fprintf(in,"  C                                With constraints : ");
	if (!opt->constraint) fprintf(in,"No\n");
	else {
         fprintf(in,"Yes\n");            
    }
    fprintf(in,"  T                        Lower bound for the rate : %f\n",opt->delta);
	fprintf(in,"  V                                  With variances : ");
	if (!opt->variance) fprintf(in,"No\n");
	else {
         fprintf(in,"Yes\n");		
         fprintf(in,"  B                          Parameter of variances : %d\n",opt->c);
         fprintf(in,"  S                                 Sequence Length : %i\n",opt->seqLength);
    }
       fprintf(in,"  R                               Estimate the root : ");              
       if (opt->estimate_root==NULL){
              fprintf(in,"No\n");                                          
       }    
       else if (strcmp(opt->estimate_root,"l")==0){
              fprintf(in,"Around the given root\n");                                          
       }           
       else if (strcmp(opt->estimate_root,"a")==0 && opt->constraint){  
            fprintf(in,"Use fast method to search on all branches\n");      
       } 
       else if (strcmp(opt->estimate_root,"a")==0 && !opt->constraint){  
            fprintf(in,"Search on all branches\n");      
       } 
       else if (strcmp(opt->estimate_root,"as")==0){  
            fprintf(in,"Use constrained mode on all branches\n");      
       } 
    fprintf(in,"  W                         Given substitution rate : ");
    if (opt->rate==NULL) fprintf(in,"No\n");
    else fprintf(in,"%s\n",opt->rate);
    if (opt->fnOutgroup==NULL)
       fprintf(in,"  G                                Remove outgroups : No\n");
    else {
         fprintf(in,"  G                         File contains outgroups : %s\n",opt->fnOutgroup);
         //opt->inFileOp = getDefaultInFileOpName(opt->inFile);
    }
	fprintf(in,"  N                               Multiple data set : ");
	if( opt->nbData< 2 )
		fprintf(in,"No\n");
	else
		fprintf(in,"Yes, %i data sets\n",opt->nbData);
    fprintf(in,"\n  H to print Help ");
	fprintf(in,"\n  Y to accept or type a letter to change an option (Q = Exit) ");
}

void setDefaultOptions( options* opt)
{
	opt->inFile = NULL;
	opt->inFileOp = NULL;
	opt->inDateFile = NULL;
	opt->outFile = NULL;
	opt->treeFile1 = NULL;	
	opt->treeFile2 = NULL;	
	opt->treeFile3 = NULL;	    	
	opt->fnOutgroup = NULL; 
	opt->seqLength = 1000;      
	opt->nbData = 1;      
    opt->rate = NULL; 
    opt->relative = false;   
    opt->mrca=0;
    opt->leaves=1;
	opt->estimate_root = NULL;
	opt->constraint = false;
    opt->variance = false;	
	opt->c = 10;
	opt->delta = 0.00001;
}

void printHelp( void )
{
	printf(BOLD"LSD: LEAST-SQUARES METHODS TO ESTIMATE RATES AND DATES - " VERSION" by Thu-Hien To\n\n");
	printf(BOLD"DESCRIPTION\n"
	FLAT"\tThis program estimates the rate and the internal dates of the input phylogenies with dated tips.\n"
    FLAT"\tIt minimizes the square errors of the branch lengths under normal distribution model.\n\n"
	);	
	printf(BOLD"SYNOPSIS\n"
	FLAT"\t" BOLD"./lsd " FLAT"[" BOLD"-i " LINE"inputFile" FLAT"] "
	FLAT"[" BOLD"-d " LINE"inputDateFile" FLAT"] "
	FLAT"[" BOLD"-o " LINE"outputFile" FLAT"] "
	FLAT"[" BOLD"-c" FLAT"] "
	FLAT"[" BOLD"-v" FLAT"] "
	FLAT"[" BOLD"-s " LINE"sequenceLength" FLAT"] "
	FLAT"[" BOLD"-n " LINE"datasetNumber" FLAT"]\n"
	FLAT"\t     [" BOLD"-t " LINE"lowerBoundRate" FLAT"] "
	FLAT"[" BOLD"-r " LINE"rootingMethod" FLAT"] "
	FLAT"[" BOLD"-b " LINE"varianceParameter" FLAT"] "
	FLAT"[" BOLD"-w " LINE"givenRateFile" FLAT"] "
	FLAT"[" BOLD"-g " LINE"outgroupFile" FLAT"] "
	FLAT"[" BOLD"-h" FLAT"]\n"
	FLAT"\n");
	
	printf(BOLD"OPTIONS\n"
	FLAT"\t" BOLD"-a " LINE"root date\n"
	FLAT"\t   If the dates of all tips are equal (which is given by option -z), you must use this option to provide the root date.\n"
	FLAT"\t   In this case, the input date file can be omitted, and the program estimates only the relative dates based on the given\n" 
    FLAT"\t   root date and tips date. By default, T[root]=0 and T[tips]=1.\n"      
    FLAT"\t" BOLD"-b " LINE"parameter of variances\n"
    FLAT"\t   The parameter to compute the variances in option -v. By default b=10.\n"
    FLAT"\t" BOLD"-c " LINE"constraints\n"    
    FLAT"\t   By using this option, we impose the constraints that the date of every node is equal or smaller then\n"
    FLAT"\t   the dates of its descendants. Without constraints, the runtime is linear (LD). With constraints, the\n"
    FLAT"\t   problem is a quadratic programming and is solved efficiently by the active-set method.\n" 
	FLAT"\t" BOLD"-d " LINE"inputDateFile\n"
	FLAT"\t   By using this options, the program reads the name of the input date file, and the file should have the following format:\n"
	FLAT"\t      n\n"
	FLAT"\t      TAXON1    DATE1\n"
	FLAT"\t      TAXON2    DATE2\n"
	FLAT"\t      ...\n"
	FLAT"\t      TAXONn    DATEn\n"
	FLAT"\t  If this option is omitted, the program will estimate relative dates by giving T[root]=0 and T[tips]=1.\n"
	FLAT"\t" BOLD"-g " LINE"outgroupFile\n"
	FLAT"\t   If your data contain outgroups, specify the name of the outgroup file here.\n" 
    FLAT"\t   The program will remove the outgroups from the trees and take the ingroup trees as the input.\n" 
    FLAT"\t   The format of this file should be:\n"
    FLAT"\t        n\n"
    FLAT"\t        OUTGROUP1\n"
    FLAT"\t        OUTGROUP2\n"
    FLAT"\t        ...\n"
    FLAT"\t        OUTGROUPn\n"     
	FLAT"\t" BOLD"-h " LINE"help\n" 
	FLAT"\t   Print this message.\n"
	FLAT"\t" BOLD"-i " LINE"inputTreesFile\n"
	FLAT"\t   The name of the input trees file. It contains tree(s) in newick format. Note that the taxa sets of all trees must be the same.\n"
	FLAT"\t" BOLD"-n " LINE"datasetNumber\n"
	FLAT"\t   The number of trees that you want to read.\n"
    FLAT"\t" BOLD"-o " LINE"outputFile\n"
	FLAT"\t   The name of the output file to write the results.\n"
    FLAT"\t" BOLD"-r " LINE"rootingMethod\n"
    FLAT"\t   This option is used to specify the rooting method to estimate the position of the root for unrooted trees, or\n"
    FLAT"\t   re-estimate the root for rooted trees. The principle is to search for the position of the root that minimizes\n" 
    FLAT"\t   the objective function.\n"
    FLAT"\t   If the tree is rooted, then either using operand \"l\" for searching the root around the given root, or using \"a\" for\n"
    FLAT"\t   searching the root on all branches. Moreover, when the constrained mode is chosen (option -c), method \"a\" firstly\n"
    FLAT"\t   estimates the root without using the constraints. After that, it uses the constrained mode to improve locally the position\n"
    FLAT"\t   of the root around this pre-estimated root. To use constrained mode on all branches in this case, please specify \"as\".\n"
    FLAT"\t   If the tree is not rooted, then the program searches the root on all branches. Similarly for the previous case, if\n"
    FLAT"\t   the constrained mode is chosen, method \"a\" uses only constrained mode to improve the root position around the pre-estimated\n"
    FLAT"\t   root which is computed without constraints. To use constraints on all branches, use \"as\".\n"
	FLAT"\t" BOLD"-s " LINE"sequenceLength\n"
	FLAT"\t   This option is used to specify the sequence length to compute the variances in the option -v.\n" 
    FLAT"\t   By default it is 1000.\n"
    FLAT"\t" BOLD"-t " LINE"lower bound for the rate\n"
    FLAT"\t   This option corresponds to the lower bound for the estimating rate.\n"
    FLAT"\t   It is 0.00001 by default.\n"
    FLAT"\t" BOLD"-v " LINE"variances\n"
    FLAT"\t   Use this option if you want to apply variances for the branch lengths in order to recompense big errors \n"
    FLAT"\t   on long estimated branch lengths. The variance of the branch Bi is Vi = (Bi+b/seq)/seq where b is specified\n"
    FLAT"\t   by option -b and seq is specified by option -s\n"
    FLAT"\t" BOLD"-w " LINE"given rate\n"
    FLAT"\t   This option is used to specify the name of the file containing the substitution rates.\n" 
    FLAT"\t   In this case, the program will use the given rates to estimate the dates of the internal nodes.\n"
    FLAT"\t   This file should have the following format\n" 
    FLAT"\t        RATE1\n"
    FLAT"\t        RATE2\n"
    FLAT"\t        ...\n"      
    FLAT"\t  where RATEi is the rate of the tree i in the inputTreesFile.\n"      
	FLAT"\t" BOLD"-z " LINE"tips date\n"
	FLAT"\t   This option is used to give the date of the tips when they are all equal. It must be used with option -a to give the\n"	
	FLAT"\t   root date. In this case the input date file can be omitted, and the program estimates only the relative dates based on\n"
	FLAT"\t   the given root date and tips date. By default, T[root]=0 and T[tips]=1.\n"    
	);
}

char* getInputString(string msg)
{
	char* outfile = (char *) malloc( sizeof(char)*MAX_FILE_NAME );
    cout<<msg<<endl;
	fgets( outfile, MAX_FILE_NAME, stdin );
	chomp( outfile );
	return outfile;
}

char* getInputTreeFileName( string msg , bool & rooted)
{
	char* outfile;
	do
	{
		outfile = getInputString(msg);
		if( access(outfile, F_OK)==0 )
			break;
		printf( "The file \"%s\" does not exist.\n", outfile );
	} while( true );
	if( access(outfile, R_OK)!=0 )
		myExit("Could not access to the file named \"%s\" in reading.\n", outfile );
    int n,m;
    counting(outfile,n,m,1);
    if (m!=2*n && m!=2*n+1) {
        cout<<"Errors in the input trees"<<endl;
        exit(EXIT_FAILURE);
    }
    rooted=true;//rooted=true if input trees are rooted
    if (m==2*n+1) rooted=false;
	return outfile;
}

char* getInputFileName( string msg)
{
	char* outfile;
	do	
	{
		outfile = getInputString(msg);
		if( access(outfile, F_OK)==0 || string(outfile).length()==0)
			break;
		printf( "The file \"%s\" does not exist.\n", outfile );
	} while( true );
	if( access(outfile, R_OK)!=0 && string(outfile).length()!=0)
		myExit("Could not access to the file named \"%s\" in reading.\n", outfile );
	return outfile;
}

char* getOutgroupFileName( string msg)
{
	char* outfile;
	do
	{
		outfile = getInputString(msg);
		if( access(outfile, F_OK)==0 || string(outfile).compare("")==0)
			break;
		printf( "The file \"%s\" does not exist.\n", outfile );
	} while( true );
    if (string(outfile).compare("")==0) {
        return outfile;
    }
	if( access(outfile, R_OK)!=0 )
		myExit("Could not access to the file named \"%s\" in reading.\n", outfile );
	return outfile;
}


void chomp( char* in )
{
	in += strlen( in )-1;
	if( *in=='\n' )
		*in='\0';
}

double getInputReal( string msg )
{
	char* word;
	do
	{
		word = getInputString( msg );
		if( isReal(word) )
			break;
		myErrorMsg("Your word is not recognized as a real.\n");
	} while( true );
	return atof( word );
}

int getInputInteger( string msg )
{
	char* word;
	do
	{
		word = getInputString( msg );
		if( isInteger(word) )
			break;
		myErrorMsg("Your word is not recognized as an integer.\n");
	} while( true );
	return atoi( word );
}

int getPositiveInputInteger( string msg )
{
	int i;
	do
	{
		i = getInputInteger(msg);
		if( i>0 )
			break;
		myErrorMsg("It must be a strictly positive integer.\n");
	} while( true );
	return i;
}

bool isOptionActivate( options* opt, char l )
{
	switch(l)
	{
		case 'i':
		case 'I':
		case 'd':
		case 'D':
		case 'o':
		case 'O':             
        case 's':
        case 'S':
		case 'c':
		case 'C':
        case 'v':
        case 'V':     
		case 'b':
		case 'B':
		case 'r':
		case 'R':
		case 'g':
		case 'G':
        case 't':
        case 'T':
        case 'w':
        case 'W':
		case 'n':
		case 'N':
		case 'y':
		case 'Y':
		case 'q':
		case 'Q':
        case 'h':
        case 'H':
			return true;
	}
	return false;
}

void setOptionsWithLetter( options* opt, char letter, bool & inputrooted, bool & out )
{
	char* fnOut;
	switch( letter )
	{
		case 'q':
		case 'Q':
			exit( EXIT_SUCCESS );
		case 'i':
		case 'I':
			free( opt->inFile );
			opt->inFile = getInputTreeFileName("Enter your Input File name> ", inputrooted);
			break;		
		case 'd':
		case 'D':
             free( opt->inDateFile );
             opt->inDateFile = getInputFileName("Enter your Input Date File name if there is any, otherwise press Enter> ");
             cout<<endl;
             if (string(opt->inDateFile).length()==0){               
                cout<<"There is no date file, so the program will estimate relative dates with root date = 0 and tips date = 1."<<endl;
	            cout<<"Type 'y' to continue or 'n' to modify the root date and the tips date"<<endl;		
	            char let[3];
	            do {
		           fgets( let, 3, stdin );	   
		           if (*let=='n' || *let=='N'){
			       do {           
				      opt->mrca = getInputReal("Enter the root date (default=0)> ");
				      opt->leaves = getInputReal("Enter the tips date (default=1)> ");        
				      if (opt->leaves <= opt->mrca) cout<<"Root date must be smaller than the tips date."<<endl;
			       } while (opt->leaves <= opt->mrca);
		          }
		          else if (*let=='y' || *let=='Y'){
				       opt->mrca=0;
				       opt->leaves=1;
		          }
		          else {
			           cout<<"Type 'y' to continue or 'n' to modify the root date and the tips date"<<endl;		
                  }		
              } while (*let!='n' && *let!='N' && *let!='y' && *let!='Y');
	          opt->relative=true;
           }
           else opt->relative=false;
	       break;
		case 's':
		case 'S':
			opt->seqLength = getPositiveInputInteger("Enter your sequence length> ");
			break;
		case 'n':
		case 'N':
			opt->nbData = getPositiveInputInteger("Enter your number of dataset> ");
			break;
		case 'o':
		case 'O':
			fnOut=getInputString("Enter your output file name > ");
            while( access( fnOut, F_OK )==0){
                cout<<"File "<<fnOut<<" already exists. Do you want to overwrite it? Y/N"<<endl;
                char letter[3];
                fgets( letter, 3, stdin );
                if (*letter=='N' || *letter=='n'){
                    fnOut = getInputString("Enter your output file name > ");
                }
                if (*letter=='Y' || *letter=='y') {
                    break;
                }
            }
            opt->outFile = new char[strlen(fnOut)];
            strcpy(opt->outFile,fnOut);
            delete[] fnOut;
            break;
		case 'c':
        case 'C':
            opt->constraint=!opt->constraint;
            break;
        case 'v':
        case 'V':
             opt->variance=!opt->variance; 
             break;
        case 'b':
        case 'B':
             if (opt->variance) opt->c = getPositiveInputInteger("Enter the parameter for the variances> ");
             break;
        case 'r':
        case 'R':
             if (inputrooted || out){ 
                if (opt->estimate_root==NULL){
                   opt->estimate_root=new char[2];
                   strcpy(opt->estimate_root,"l");
                }
                else if (strcmp(opt->estimate_root,"l")==0){
                   strcpy(opt->estimate_root,"a");
                }
                else if (strcmp(opt->estimate_root,"a")==0 && opt->constraint){
                   strcpy(opt->estimate_root,"as");
                }
                else if (strcmp(opt->estimate_root,"a")==0 && !opt->constraint) opt->estimate_root=NULL;
                else if (strcmp(opt->estimate_root,"as")==0) opt->estimate_root=NULL;
             }
             else {
                  if (opt->estimate_root==NULL){
                     opt->estimate_root=new char[2];
                     strcpy(opt->estimate_root,"a");
                  }   
                  else if (strcmp(opt->estimate_root,"a")==0 && opt->constraint){
                     strcpy(opt->estimate_root,"as");
                  }
                  else if (strcmp(opt->estimate_root,"a")==0 && !opt->constraint){
                     cout<<"The trees are not rooted, you must use either option -g to specify the outgroups file or -r to estimate the root"<<endl;
                  }
                  else if (strcmp(opt->estimate_root,"as")==0){
                     strcpy(opt->estimate_root,"a");
                  }
             }
             break;
		case 'g':
		case 'G':
            if (opt->fnOutgroup==NULL){
                char* fo = getOutgroupFileName("Enter the name of the file that contains your outgroups> ");
                if (string(fo).compare("")!=0) {
                    out=true;
                    opt->fnOutgroup=fo;
                }
            }
            else{
                opt->fnOutgroup=NULL;
                out=false;
                if (!inputrooted && opt->estimate_root==NULL){
                   opt->estimate_root=new char[2];
                   strcpy(opt->estimate_root,"a");              
                } 
            }
            break;
		case 'w':
        case 'W':
             if (opt->rate==NULL) {
                opt->rate = getInputFileName("Enter the name of the file that contains the rates> ");				              
                if (string(opt->rate).length()==0) opt->rate=NULL;
             }
             else opt->rate=NULL;
             break;
		case 't':
        case 'T':
             opt->delta = getInputReal("Enter the lower bound for the rate> ");
             break;
        case 'h':
        case 'H':
            printHelp();
            break;
	}
}


