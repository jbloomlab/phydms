#ifndef OPTIONS_H
#define OPTIONS_H
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "stdarg.h"
#include <math.h>
#include "utils.h"
//Pour affiche l'aide
#ifndef VERSION
#define VERSION "v0.2"
#endif
#ifndef BOLD
#define BOLD      "\033[00;01m"
#endif
#ifndef FLAT
#define FLAT      "\033[00;00m"
#endif
#ifndef LINE
#define LINE      "\033[00;04m"
#endif
//Longuer Max d'un nom de fichier
#ifndef MAX_FILE_NAME
#define MAX_FILE_NAME 256
#endif
#ifndef MAX_NAME_LENGTH
#define MAX_NAME_LENGTH 256
#endif


using namespace std;

typedef struct options
{
	char*  inFile;         //Nom du fichier d'arbres
	char* inFileOp;        //Nom du fichier d'arbres apres enlevant outgroup
	char*  inDateFile;     //Nom du fichier contenant les dates
	char*  outFile;        //Nom du fichier de resultats.
	char* treeFile1;       //Nom du fichier d'abres sorties Nexus
	char* treeFile2;       //Nom du fichier d'arbres sorties Newick avec des longueurs de branches mesures par substiution par site
	char* treeFile3;       //Nom du fichier d'arbres sorties Newick avec des longueurs de branches mesures par annees
	bool relative;         //=true if all the leaves have the same date, the program estimate the relative dates
	double mrca;
	double leaves;
	int    seqLength;      //Longueur des sequences dans l'alignement
	int    nbData;         //Nombre de cas a  traiter (dans le cas de bootstrap)
	char* fnOutgroup;
	char* rate;           //le fichier contient les taux en entree
//	string estimate_root;
	char* estimate_root;    //Method to estimate root 
	bool constraint;       //Impose the constraints or not
	bool variance;         //Use the variances or not
	int  c;                //var = b+c/s;	
	double delta;             //the lower bound of the rate
} options;


options* getOptions( int, char** , bool&);
options* getCommandLine( int, char**, bool&);
options* getInterface( bool & );
void     printHelp( void );
void     setDefaultOptions( options* );
char*    getInputTreeFileName( string , bool&);
char*    getInputFileName( string );
char*    getOutgroupFileName( string );
void     chomp( char* );
void     printInterface( FILE*, options*, bool);
char*    getDefaultOutputFileName( string );
char*    getDefaultOutputNexusTreeFileName( string );
char*    getDefaultOutputNewick1TreeFileName( string );
char*    getDefaultOutputNewick2TreeFileName( string );
char*    getDefaultInFileOpName( string );
void     setOptionsWithLetter( options* , char , bool &, bool &);
double   getInputReal( string );
int      getInputInteger( string );
int      getPositiveInputInteger( string );
char*    getInputString( string );
bool     isOptionActivate( options*, char );
FILE*    openOutputFile( char** );
FILE*    myFopen( char*, char* );
#endif

