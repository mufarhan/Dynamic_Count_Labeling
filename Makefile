CC = g++ -std=c++2a -O3 -Wall -Wextra -pthread -o
TCC = g++ -std=c++2a -ggdb -Wall -Wextra -o
INC = src/road_network.cpp src/util.cpp

all: index query update

index:
	$(CC) index src/index.cpp $(INC)
query:
	$(CC) query src/query.cpp $(INC)
update:
	$(CC) update src/update.cpp $(INC)

clean:
	rm index query update

.PHONY: index query update
