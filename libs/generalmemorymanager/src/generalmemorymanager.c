#include "generalmemorymanager.h"
#include <stdlib.h>
#include <stdio.h>

//Debug?:
//#define DEBUG_CTRL
//macros:
/**
 * @def ERROR(message)
 *
 * Displays "message" as error-message with
 * general info.
 *
 */
#define ERROR(message) fprintf(stderr,"%s : %s in %s (Line: %d)\n",\
						(message),__func__,__FILE__,__LINE__)

/**
 * @def MAXTRIES
 *
 * How many Tries will be attempted for an action to take place
 *
 */
#define MAXTRIES 5

void* mem_init(int showinfo){
	return mem_backbone (MEM_INITIALIZE,NULL,showinfo);
}

void* mem_alloc(size_t size){
	return mem_backbone (MEM_ALLOCATION,NULL,size);
}

void* mem_realloc(void* pointer,size_t size){
	return mem_backbone (MEM_REALLOCATION,pointer,size);
}

void mem_free(void* pointer){
	mem_backbone (MEM_FREE,pointer,0);
	return ;
}

void mem_free_all(void){
	mem_backbone (MEM_FREEALL,NULL,0);
	return;
}

void* mem_backbone(const mem_instruction_t instruction,void* pointer2process,
					size_t parameter){
	/**
	 * Declarations:
	 *
	 * @var memory_len		length of array all_mem_blocks
	 * @var amount_elements	Actual amount of elements in all_mem_blocks
	 *						(Not "NULL")
	 * @var all_mem_blocks	Pointer to array of all allocated memory blocks
	 * @var exitcode		Saves the exitcode of lower level functions
	 *						to pass it onto the top level
	 * @var try		 		Counts how many tries the action took
	 * @var retpointer		Pointer to be returned by mem_backbone
	 *
	 */
	static int memory_len;
	static int amount_elements;
	static void** all_mem_blocks;
	static int showinfo;
	/**
	 * @note Test rather this is the initial call
	 */
	if(all_mem_blocks==NULL&&instruction!=MEM_INITIALIZE){
		printf("[WARNING] Use of General memory manager without"
				" Initialization\n");
		mem_init(1);
	}
	mem_exitcode_t exitcode;
	int try=0;
	void * retpointer=NULL;

	/**
	 * @note	Tries "MAXTRIES" amount of times to fulfill the required action
	 * 			Each time the required action "instruction" is witched trough 
	 * 			and the needed lower level function is called.
	 * 			Afterwards it gets checked rather the call was successful.
	 *
	 */
	do{
		exitcode=MEM_NOPROBLEMS;
		try++;
		switch (instruction)
		{
		case MEM_INITIALIZE:
			if(parameter){
				printf("[INFO] Beginning memory manager Initialization\n");
			}
			mem_backinit(&memory_len,&amount_elements,
							&all_mem_blocks,&exitcode,&showinfo,parameter);
			retpointer=all_mem_blocks;
			break;
		case MEM_ALLOCATION:
			if(showinfo){
				printf("[INFO] Beginning memory allocation of size: %ld\n",
					parameter);
			}
			retpointer=mem_backalloc (&memory_len, &amount_elements,
							  &all_mem_blocks, &exitcode, parameter);
			break;
		case MEM_REALLOCATION:
			if(showinfo){
				printf("[INFO] Beginning memory reallocation of size %ld at"
					"%p\n", parameter,pointer2process);
			}
			retpointer=mem_backrealloc(&memory_len, &all_mem_blocks,
								&exitcode, parameter ,pointer2process);
			break;
		case MEM_FREE:
			if(showinfo){
				printf("[INFO] Beginning to free memory block at %p\n",
					pointer2process);
			}
			mem_backfree(&memory_len, &amount_elements, &all_mem_blocks, 
						 &exitcode,pointer2process);
			break;
		case MEM_FREEALL:
			if(showinfo){
				printf("[INFO] Beginning to free all memory blocks\n");
			}
			mem_backfreeall(&memory_len, &amount_elements, &all_mem_blocks, 
						 &exitcode,&showinfo);
			break;
		}
#ifdef DEBUG_CTRL
		printf("Status: \nmemory_len=%d\namount_elements"
			"=%d\nexitcode=%d\nall_mem_blocks=%p\nretpointer=%p\n",
			memory_len,amount_elements,exitcode,(void*)all_mem_blocks,
			retpointer);
#endif //ifdef DEBUG_CTRL
	}while((exitcode!=MEM_NOPROBLEMS)&&(try<MAXTRIES));
	if(try>=MAXTRIES){
		ERROR("Could not resolve previous Error!");
		exit(exitcode);
	}



	return retpointer;
}

void mem_backinit(int* memory_len_ptr,int* amount_elements_ptr,
				void*** all_mem_blocks_ptr,mem_exitcode_t* exitcode_ptr,
				int* showinfo_ptr,int parameter){
	/**
	 * @note	"memory_len,amount_elements" get initialized and 
	 *			memory gets allocate for "all_mem_blocks".
	 *			Afterwards the memory allocation for all_mem_blocks
	 *			gets registered in all_mem_blocks(allways at position 0)
	 * 
	 */
	*showinfo_ptr=1;
	*memory_len_ptr=1;
	*amount_elements_ptr=1;
	*all_mem_blocks_ptr=malloc(sizeof (void*)*(*memory_len_ptr));
	if(*all_mem_blocks_ptr==NULL){
		ERROR("Memory allocation error for all_mem_blocks");
		*exitcode_ptr=ERR_ALLOCATION;
		return ;
	}
	(*all_mem_blocks_ptr)[0]=*all_mem_blocks_ptr;
	if(parameter==0){
		*showinfo_ptr=0;
	}
	return ;
}

void* mem_backalloc(int* memory_len_ptr,int* amount_elements_ptr,
					void*** all_mem_blocks_ptr,mem_exitcode_t* exitcode_ptr,
					size_t size){
	/**
	 * Declarations:
	 *
	 * @var new_mem_block	Pointer to the newly allocated memory block
	 * @var index			Index of the new pointer in all_mem_blocks
	 * @var temp			Stores returned pointers temporarily
	 */
	void* new_mem_block=NULL;
	int index;
	void* temp;

	/**
	 * @note	check rather all_mem_blocks has an empty spot to save 
	 * 			the new_mem_block pointer and increase its size if not.
	 * 			Adjust memory_len to match current size of all_mem_blocks.
	 * 
	 */
	if(*amount_elements_ptr<*memory_len_ptr){
		index=mem_searcher(memory_len_ptr,all_mem_blocks_ptr,exitcode_ptr,
					 NULL);
	}else{
		index=*memory_len_ptr;
		temp=realloc (*all_mem_blocks_ptr,
				  (*memory_len_ptr+1)*sizeof(void*));
		if(temp==NULL){
			ERROR("Memory reallocation error for all_mem_blocks");
			*exitcode_ptr=ERR_REALLOCATION;
		}else{
			*all_mem_blocks_ptr=temp;
			(*memory_len_ptr)++;
			(*all_mem_blocks_ptr)[index]=NULL;
		}
	}
	/**
	 *@note	Allocate new memory block and save its address in 
	 * 		all_mem_blocks.
	 * 		Adjust amount_elements to match current amount of 
	 * 		elements in all_mem_blocks.
	 * 
	 */
	new_mem_block=malloc (size);
	if(new_mem_block==NULL){
		ERROR("Memory allocation error for new_mem_block");
		*exitcode_ptr=ERR_ALLOCATION;
	}
	if(*exitcode_ptr!=MEM_NOPROBLEMS){
		return NULL;
	}
	(*all_mem_blocks_ptr)[index]=new_mem_block;
	(*amount_elements_ptr)++;
	return new_mem_block;
}

void* mem_backrealloc(int* memory_len_ptr, void*** all_mem_blocks_ptr,
						mem_exitcode_t* exitcode_ptr, size_t newsize,
						void* pointer2realloc){
	/**
	 * Declarations:
	 *
	 * @var adjusted_mem_block	Pointer to the new reallocated memory block
	 * @var index		Index of "pointer2realloc" in all_mem_blocks
	 */
	void* adjusted_mem_block=NULL;
	int index;
	
	/**
	 * @note	Detect rather "pointer2realloc" is known to the 
	 * 			memorymanager and reallocated the memory block if it is
	 * 			and save the new(might be changed) address in all_mem_blocks.
	 * 
	 */
	index=mem_searcher(memory_len_ptr,all_mem_blocks_ptr,exitcode_ptr,
						pointer2realloc);
	if(*exitcode_ptr!=MEM_NOPROBLEMS){
		return NULL;
	}
	
	adjusted_mem_block=realloc (pointer2realloc,newsize);
	if(adjusted_mem_block==NULL){
		ERROR("Memory reallocation error for adjusted_mem_block");
		*exitcode_ptr=ERR_REALLOCATION;
	}else{
		(*all_mem_blocks_ptr)[index]=adjusted_mem_block;
	}
	return adjusted_mem_block;
}

void mem_backfree(int* memory_len_ptr,int* amount_elements_ptr,
				  void*** all_mem_blocks_ptr,mem_exitcode_t* exitcode_ptr,
				  void* pointer2free){
	/**
	 * Declarations:
	 *
	 * @var index		Index of "pointer2free" in all_mem_blocks
	 */
	int index;
	
	/**
	 * @note	Search for "pointer2free" and free its memory block,
	 * 			afterwards adjust "amount_elements_ptr" to match the 
	 * 			current amount of elements saved in "all_mem_block".
	 */
	index=mem_searcher(memory_len_ptr,all_mem_blocks_ptr,exitcode_ptr,
						pointer2free);
	if(*exitcode_ptr!=MEM_NOPROBLEMS){
		return;
	}
	free(pointer2free);
	(*all_mem_blocks_ptr)[index]=NULL;
	(*amount_elements_ptr)--;
	return;
}

void mem_backfreeall(int* memory_len_ptr,int* amount_elements_ptr,
				  void*** all_mem_blocks_ptr,mem_exitcode_t* exitcode_ptr,
				  int* showinfo_ptr){
	/**
	 * Declarations:
	 *
	 * @var i		Iteration variable (returned index)
	 *
	 */
	int i;
	/**
	 * @note	Free all memory blocks (except all_mem_blocks itself) and check 
	 * 			rather there were to many or to less pointers.
	 * 			If everything was fine also free all_mem_blocks itself.
	 */
	for(i=*memory_len_ptr-1;i>0;i--){
		if((*all_mem_blocks_ptr)[i]==NULL){
			continue;
		}
		free((*all_mem_blocks_ptr)[i]);
		(*all_mem_blocks_ptr)[i]=NULL;
		(*amount_elements_ptr)--;
	}
	if((*amount_elements_ptr)>1){
		ERROR("Freed less pointers than expected");
		*exitcode_ptr=ERR_POINTERMISSING;
		return;
	}
	if((*amount_elements_ptr)<1){
		ERROR("Freed more pointers than expected");
		*exitcode_ptr=ERR_POINTEROVERHEAD;
		return;
	}
	free(*all_mem_blocks_ptr);
	*all_mem_blocks_ptr=NULL;
	*amount_elements_ptr=0;
	*memory_len_ptr=0;
	if(*showinfo_ptr){
		printf("[INFO] Freed all memory blocks\n");
	}
	return;
}

int mem_searcher(int* memory_len_ptr, void*** all_mem_blocks_ptr,
				mem_exitcode_t* exitcode_ptr, void* pointer2search){
	/**
	 * Declarations:
	 *
	 * @var i	Iteration variable (returned index)
	 *
	 */
	int i;
	
	/**
	 * @note Search for pointer2search, if it is not found 
	 * "ERR_POINTERNOTFOUND" is raised.
	 */
	for(i=*memory_len_ptr-1;i>=0;i--){
		if((*all_mem_blocks_ptr)[i]==pointer2search){
			break;
		}
	}
	if(i<0){
		ERROR("Could not find specified pointer");
		*exitcode_ptr=ERR_POINTERNOTFOUND;
	}
	return i;
}
#undef ERROR
#undef DEBUG_CTRL
