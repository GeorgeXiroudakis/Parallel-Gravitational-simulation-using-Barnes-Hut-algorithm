all: ask3

ask3: src/ask3.cpp
	g++ src/ask3.cpp -o ask3 -ltbb

clean:
	@ rm -f ask3 output.txt
