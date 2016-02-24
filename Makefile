CC=g++
BOOST_LIB=/usr/lib64

fixedTime : driveFixedTime.o experiment.o population.o organism.o parameters.o rv_generators.o
	${CC} -o fixedTime.exe  driveFixedTime.o experiment.o population.o organism.o parameters.o rv_generators.o -L ${BOOST_LIB} -lboost_serialization -L /usr/include -lgsl ${MKL} -Wall -O3

compete : driveCompete.o experiment.o population.o organism.o parameters.o rv_generators.o
	${CC} -o compete.exe  driveCompete.o experiment.o population.o organism.o parameters.o rv_generators.o -L ${BOOST_LIB} -lboost_serialization -L /usr/include -lgsl ${MKL} -Wall -O3
	
printCompete : driveCompetePrint.o experiment.o population.o organism.o parameters.o rv_generators.o
	${CC} -o printCompete.exe  driveCompetePrint.o experiment.o population.o organism.o parameters.o rv_generators.o -L ${BOOST_LIB} -lboost_serialization -L /usr/include -lgsl ${MKL} -Wall -O3
	                      
driveFixedTime.o: driveFixedTime.cpp experiment.o population.o organism.o parameters.o rv_generators.o 
	${CC} -c -I${BOOST_LIB} driveFixedTime.cpp -Wall -O3 -o driveFixedTime.o

driveCompete.o: driveCompete.cpp experiment.o population.o organism.o parameters.o rv_generators.o 
	${CC} -c -I${BOOST_LIB} driveCompete.cpp -Wall -O3 -o driveCompete.o
	
driveCompetePrint.o: driveCompetePrint.cpp experiment.o population.o organism.o parameters.o rv_generators.o 
	${CC} -c -I${BOOST_LIB} driveCompetePrint.cpp -Wall -O3 -o driveCompetePrint.o
	
experiment.o: experiment.cpp experiment.hpp population.o rv_generators.o 
	${CC} -c experiment.cpp -I${BOOST_LIB} -O3 -Wall

population.o: population.cpp population.hpp temp_templates.hpp organism.o rv_generators.o
	${CC} -c population.cpp -I${BOOST_LIB} -O3 -Wall

organism.o: organism.cpp organism.hpp temp_templates.hpp parameters.o 
	${CC} -c organism.cpp -I${BOOST_LIB} -O3 -Wall

parameters.o: parameters.cpp parameters.hpp temp_templates.hpp 
	${CC} -c parameters.cpp -I${BOOST_LIB} -O3 -Wall

rv_generators.o: rv_generators.cpp rv_generators.hpp
	${CC} -c rv_generators.cpp -I${BOOST_LIB} -O3 -Wall
	
clean:
	rm *.o
