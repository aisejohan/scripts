test:
	gcc -Wall -O3 -march=native -c scalar.c compute.c pol.c helper.c test_scalars.c grobner.c
	gcc -Wall -O3 -march=native -o tester test_scalars.o compute.o pol.o helper.o scalar.o grobner.o

input_pol:
	gcc -Wall -DINPUT_F -O3 -march=native -c scalar.c compute.c pol.c helper.c test_scalars.c grobner.c
	gcc -Wall -DINPUT_F -O3 -march=native -o tester test_scalars.o compute.o pol.o helper.o scalar.o grobner.o

clean:
	rm -f tester
	rm -f *.o

debug:
	gcc -g -Wall -pedantic -std=c99 -march=native -c scalar.c compute.c pol.c helper.c test_scalars.c grobner.c
	gcc -g -Wall -march=native -o tester test_scalars.o compute.o pol.o helper.o scalar.o grobner.o

debug_input:
	gcc -g -Wall -DINPUT_F -pedantic -std=c99 -march=native -c scalar.c compute.c pol.c helper.c test_scalars.c grobner.c
	gcc -g -Wall -DINPUT_F -march=native -o tester test_scalars.o compute.o pol.o helper.o scalar.o grobner.o

output_pol:
	gcc -DINPUT_F -DOUTPUT_LIST -Wall -c scalar.c compute.c pol.c helper.c test_scalars.c grobner.c
	gcc -Wall -o tester scalar.o compute.o pol.o helper.o test_scalars.o grobner.o
