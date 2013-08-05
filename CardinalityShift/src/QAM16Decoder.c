

unsigned int decodeQAM16(float* r, float* dist)
{
  int i=0;
  while(i<17 && r[0]> i* .25-1){
      i++;
  }
  int j = 0;
  while(j<17 && r[1]> j* .25-1){
      j++;
  }
  return (i-1)+16*(j-1);
}

