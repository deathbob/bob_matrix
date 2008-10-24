all:
	g++ one_417.cpp -o 417one -g -Wall -lm
test:
	g++ 417_test.cpp -o 417test -g -Wall -lm