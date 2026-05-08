CXX = clang++
CFLAGS = -I./pml/ntl-extras/ -I./pml/ntl-extras/util/ -I./pml/ntl-extras/include
LDFLAGS = ./pml/ntl-extras/lib/*.o -lntl -pthread

main.o: main.cpp pml
	$(CXX) -O3 -march=native main.cpp -o main $(CFLAGS) $(LDFLAGS)

pml:
	git clone https://github.com/vneiger/pml
	cd pml/ntl-extras && make

clean:
	rm main