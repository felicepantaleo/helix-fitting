#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
/*
void mypointer()
{
	cout << "mypointer" << endl;
};
*/

void mypointer(int a)
{
	printf("mypointer mit arg");
};

// void call_pointer ( void (*pointer)() )
void call_pointer ( void (*pointer)(int) )
{
	printf("%p %p\n",mypointer, pointer) ;
	pointer(3);
};

int main (int argc, char **argv)
{
	int h=3;
	printf("okay");
	call_pointer( mypointer );
};
