
#define CHECK_NEW_INT(size) \
	{
		
	}
#define	malloc_no_free(type, type_str)						    \
		{									    \
		printf("%s\tsizeof(%s)=%d\t", type_str, type_str, sizeof(type));\
		char i;							    \
		type *pthis = NULL, *plast = NULL;				    \
		for(i = 0; i < 6; ++i){					    \
		if( NULL != (pthis = (type *)malloc(sizeof(type))) ){   \
		if(plast)					    \
		printf("%d\t", (int)pthis - (int)plast);\
		plast	= pthis;/*a litte logic to p*/      \
		/*free(pthis);*/				    \
		pthis	= NULL;				    \
		}							    \
		}								    \
		printf("\n\n");						    \
		}