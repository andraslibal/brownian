alpha:
	g++ alpha.cpp -o alpha.out
	./alpha.out test.txt

clean:
	rm alpha.out test.txt