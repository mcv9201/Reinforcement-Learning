#include <stdio.h>
#include<math.h>
#include<stdlib.h>

float field(int x, int y)
{
    float result;
    float pi = 22/7;
    result = 0.0132 + (-0.098140)*cos(2*pi*y/100) + (0.444264)*sin(2*pi*y/100) + (-0.024930)*cos(2*pi*x/100) +
              (0.528388)*sin(2*pi*x/100) + (0.409736)*cos(2*pi*(x+y)/100) + (0.484270)*sin(2*pi*(x+y)/100);
    float new_res;
    new_res = 1.0/(1+exp(-2.0*(result - 0.2)));
    return new_res;
}


int main()
{
  FILE *fpw;
  fpw = fopen("field.txt","w");
  for(float i=-50;i<=50;i=i+0.5)
  {
    for(float j=-50;j<=50;j=j+0.5)
    {
      fprintf(fpw,"%f\t%f\t",i,j);
      fprintf(fpw,"%f\n",field(i,j));
    }
  }
  fclose(fpw);
  return 0;
}
