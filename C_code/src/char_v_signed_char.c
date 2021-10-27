#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <math.h>

// char vs. signed char on this system
// Compile with - gcc char_v_signed_char.c -o char_v_schar.exe -lm
int main(){
	char* data_char;
	signed char* data_schar;
	data_char = (char*)calloc(2, sizeof(char));
	data_schar = (signed char*)calloc(2, sizeof(signed char));	
	data_char[0] = 3.22;
	data_char[1] = -3.22;
	data_schar[0] = 3.22;
	data_schar[1] = -3.22;

	float data_float = data_char[0]*1.0f;
	float data_float2 = data_char[1]*1.0f;

	float data_float3 = data_char[0];
	float data_float4 = data_char[1];

	printf("data_char[0] =  %d and data_char[1] = %d\n", data_char[0], data_char[1]);
	printf("data_schar[0] =  %d and data_schar[1] = %d\n", data_schar[0], data_schar[1]);
	printf("data_float =  %0.1f and data_float2 =  %0.1f \n", data_float, data_float2);
	printf("data_float3 =  %0.1f and data_float4 =  %0.1f \n", data_float3, data_float4);
	
	return 0;
}
