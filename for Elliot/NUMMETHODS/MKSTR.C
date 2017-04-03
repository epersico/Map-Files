#include <stdio.h>
#include <string.h>

#define NL(flnm,dumch) {do dumch=fgetc(flnm); while (dumch!='\n' && dumch!=EOF);}


void mkstr(char *cnummer, unsigned long int nummer)
{
    unsigned long int tmp=nummer;
    char tmpstr[160], nostr[160];

    strcpy(tmpstr,"");
    while (tmp>0) {
	switch (tmp-10*(tmp/10)) {
	    case 1: {strcpy(nostr,"1"); break;}
	    case 2: {strcpy(nostr,"2"); break;}
	    case 3: {strcpy(nostr,"3"); break;}
	    case 4: {strcpy(nostr,"4"); break;}
	    case 5: {strcpy(nostr,"5"); break;}
	    case 6: {strcpy(nostr,"6"); break;}
	    case 7: {strcpy(nostr,"7"); break;}
	    case 8: {strcpy(nostr,"8"); break;}
	    case 9: {strcpy(nostr,"9"); break;}
	    case 0: {strcpy(nostr,"0"); break;}
	}
	strcat(nostr,tmpstr);
	strcpy(tmpstr,nostr);
      	tmp/=10;
    }
    strcpy(cnummer,tmpstr);
}
