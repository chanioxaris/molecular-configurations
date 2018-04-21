OBJS = structs.o main.o input_functions.o cRMSD.o clustering_init.o clustering_assignment.o clustering_update.o functions.o output_functions.o silhouette.o metric_functions.o quicksort.o
EXEC = proteins
GCCC = gcc -c
GCCO = gcc -o
LIBS = -lm -lblas -llapacke

$(EXEC): $(OBJS)
	$(GCCO) $(EXEC) $(OBJS) $(LIBS)

structs.o: structs.c
	$(GCCC) structs.c

main.o: main.c
	$(GCCC) main.c

input_functions.o: input_functions.c
	$(GCCC) input_functions.c
	
cRMSD.o: cRMSD.c
	$(GCCC) cRMSD.c
	
clustering_init.o: clustering_init.c
	$(GCCC) clustering_init.c
	
clustering_assignment.o: clustering_assignment.c
	$(GCCC) clustering_assignment.c
	
clustering_update.o: clustering_update.c
	$(GCCC) clustering_update.c
	
functions.o: functions.c
	$(GCCC) functions.c

output_functions.o: output_functions.c
	$(GCCC) output_functions.c	
	
silhouette.o: silhouette.c
	$(GCCC) silhouette.c	

metric_functions.o: metric_functions.c
	$(GCCC) metric_functions.c	

quicksort.o: quicksort.c
	$(GCCC) quicksort.c	
	
clean:
	rm -rf $(OBJS) $(EXEC)
