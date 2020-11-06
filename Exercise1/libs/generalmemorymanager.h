/**
 *  @file generalmemorymanager.h
 *  @brief Provides a general framework for memorymanagment
 *  @author Dichter nd@nocoffeetech.de
 *
 *
 *
 *	@date 18.04.2020	First Implementation of:
 *						mem_init,mem_alloc,mem_realloc,
 *						mem_free,mem_free_all,mem_backbone,
 *						mem_backinit
 *						Prototyp of:
 *						mem_backalloc
 * 	@date 19.04.2020	First Implementation of:
 * 						mem_backalloc,mem_backrealloc,
 * 						mem_backfree,mem_backfreeall,
 * 						mem_searcher
 * 						Initial testing was successful!
 * 	@date 19.04.2020	Added check in mem_backbone rather 
 * 						all_mem_blocks is set/initialcall
 *
 *
 *  @todo add function to allocate memory block and split it into sections
 *  @todo add function to reorganize all_mem_blocks
 *
 *  @test ./test_generalmemorymanager.c
 *
 *  @bug No known bugs
 *
 *  @version 1.1
 */

#ifndef GENERALMEMORYMANAGER_H
#define GENERALMEMORYMANAGER_H
#include <stdlib.h>


/**
 * @typedef typedef enum file_error_t
 * @brief Typedefinition of a enum to respresent the \
 * instructions given to mem_backbone
 */
typedef enum mem_instruction_t {
	MEM_INITIALIZE = 0,
	MEM_ALLOCATION = 1,
	MEM_REALLOCATION = 2,
	MEM_FREE = 3,
	MEM_FREEALL = 4
} mem_instruction_t;



/**
 * @typedef typedef enum file_error_t
 * @brief Typedefinition of a enum to represent the \
 * instructions given to mem_backbone
 */
typedef enum mem_exitcode_t {
	MEM_NOPROBLEMS = 50,
	ERR_ALLOCATION = 51,
	ERR_REALLOCATION = 52,
	ERR_POINTERNOTFOUND = 53,
	ERR_POINTERMISSING = 54,
	ERR_POINTEROVERHEAD = 55,
	ERR_ALL_MEM_BLOCKS_NOT_DEFINED = 56
} mem_exitcode_t;



/**
 * @fn int mem_init(int showinfo);
 * @brief Initializes the memorymanager
 *
 * Calls mem_backbone to initialize the memory manager.
 * In particular the array of
 * pointers to allocated memory blocks gets initialized.
 *
 * @param showinfo Rather info messages should be shown
 * 
 * @return Pointer to array of Pointers to allocated memory blocks
 *
 */
void* mem_init(int showinfo);



/**
 * @fn void* mem_alloc(size_t size);
 * @brief Allocates "size" sized memory block
 *
 * Calls mem_backbone,which calls mem_backalloc to allocate
 * a memory block and make it known to the memory manager.
 *
 * @param size Size of memory block to be allocated (in Byte)
 *
 * @return Pointer to the allocated memory block
 *
 */
void* mem_alloc(size_t size);



/**
 * @fn void* mem_realloc(void* pointer,size_t size);
 * @brief Reallocates memory at "pointer" to size "size"
 *
 * Calls mem_backbone,which calls mem_backrealloc to reallocate
 * the memory block at "pointer" and
 * make it known to the memory manager if the address changed.
 *
 * @param size Size of memory block to be reallocated (in Byte)
 * @param pointer Pointer to existing memory block
 *
 * @return Pointer to the reallocated memory block
 *
 */
void* mem_realloc(void* pointer,size_t size);



/**
 * @fn void mem_free(void* pointer);
 * @brief Frees memory at "pointer"
 *
 * Calls mem_backbone,which calls mem_backfree to free
 * the memory block at "pointer" and
 * make it known to the memory manager.
 *
 * @param pointer Pointer to existing memory block to be freed
 *
 */
void mem_free(void* pointer);



/**
 * @fn void mem_free_all(void);
 * @brief Frees all memory blocks
 *
 * Calls mem_backbone,which calls mem_backfree_all to free
 * all the memory blocks saved in all_mem_blocks.(including itself)
 * Afterwards mem_init needs to be called again to use the memorymanager.
 *
 */
void mem_free_all(void);



/**
 * @fn void* mem_backbone(const mem_instruction_t instruction,\
						void* pointer, size_t parameter);
 * @brief Handles all memory interactions (should never be called directly)
 *
 * Handles all memory interactions and hands them off to mem_backinit,
 * mem_backalloc, mem_backrealloc, mem_backfree, mem_backfree_all.
 *
 * @param instruction Tells, which functionality is needed
 * @param pointer2process Pointer to be processed
 * @param parameter Size of the needed memory block or other options
 *
 * @return Pointer to the (re)allocated memory block \
			/array of pointers to all memory blocks
 */
void* mem_backbone(const mem_instruction_t instruction,void* pointer2process,
			 size_t parameter);



/**
 * @fn void mem_backinit(int* memory_len_ptr,int* amount_elements_ptr,\
 * 						void*** all_mem_blocks_ptr,int* exitcode_ptr,
 * 						int* showinfo_ptr,int parameter);
 * @brief Initializes the array of pointers to allocated memory blocks\
						(should never be called directly)
 *
 * The array of
 * pointers to allocated memory blocks gets initialized.
 *
 * @param memory_len_ptr Pointer to current len of "all_mem_blocks"\
							in the memory in units of sizeof(void*)
 * @param amount_elements_ptr Amount of memory blocks currently allocated
 * @param all_mem_blocks_ptr Array of all memory block pointers
 * @param exitcode_ptr Exitcode address (needed in later error processing)
 */
void mem_backinit(int* memory_len_ptr,int* amount_elements_ptr,
				void*** all_mem_blocks_ptr,mem_exitcode_t* exitcode_ptr,
				int* showinfo_ptr,int parameter);



/**
 * @fn void* mem_backalloc(int* memory_len_ptr,int* amount_elements_ptr,\
				void*** all_mem_blocks_ptr,int* exitcode_ptr,size_t size);
 * @brief Allocates a memory block of size "size"\
				(should never be called directly)
 *
 * Allocates a memory block of size "size" and stores its address in
 * all_mem_blocks.
 *
 * @param memory_len_ptr Pointer to current len of "all_mem_blocks"\
							in the memory in units of sizeof(void*)
 * @param amount_elements_ptr Amount of memory blocks currently allocated
 * @param all_mem_blocks_ptr Array of all memory block pointers
 * @param exitcode_ptr Exitcode address (needed in later error processing)
 * @param size Size of the memory block to be allocated (in bytes)
 *
 * @return Pointer to the allocated memory block
 */
void* mem_backalloc(int* memory_len_ptr,int* amount_elements_ptr,
				void*** all_mem_blocks_ptr,mem_exitcode_t* exitcode_ptr,
				size_t size);



/**
 * @fn void* mem_backrealloc(int* memory_len_ptr, void*** all_mem_blocks_ptr,\
				mem_exitcode_t* exitcode_ptr, size_t newsize,\
				void* pointer2realloc);
 * @brief Reallocates the memory block at "pointer2realloc" to size "newsize"\
						(should never be called directly)
 *
 * Reallocates the memory block at "pointer2realloc" to size "newsize"
 * and stores its address in all_mem_blocks.
 *
 * @param memory_len_ptr Pointer to current len of "all_mem_blocks"\
							in the memory in units of sizeof(void*)
 * @param all_mem_blocks_ptr Array of all memory block pointers
 * @param exitcode_ptr Exitcode address (needed in later error processing)
 * @param newsize New size of the memory block to be reallocated (in bytes)
 * @param pointer2realloc Pointer to memory block to be reallocated
 *
 * @return Pointer to the new reallocated memory block
 */
void* mem_backrealloc(int* memory_len_ptr, void*** all_mem_blocks_ptr,
				mem_exitcode_t* exitcode_ptr, size_t newsize,
				void* pointer2realloc);



/**
 * @fn void mem_backfree(int* amount_elements_ptr,void*** all_mem_blocks_ptr,\
						mem_exitcode_t* exitcode_ptr,void* pointer2free)
 * @brief Frees the memory block at "pointer2free"\
				(should never be called directly)
 *
 * Frees the memory block at "pointer2free" 
 * and deletes it from "all_mem_blocks".
 *
 * @param memory_len_ptr Pointer to current len of "all_mem_blocks"\
							in the memory in units of sizeof(void*)
 * @param amount_elements_ptr Amount of memory blocks currently allocated
 * @param all_mem_blocks_ptr Array of all memory block pointers
 * @param exitcode_ptr Exitcode address (needed in later error processing)
 * @param pointer2realloc Pointer to memory block to be freed
 */
void mem_backfree(int* memory_len_ptr,int* amount_elements_ptr,
				  void*** all_mem_blocks_ptr,mem_exitcode_t* exitcode_ptr,
				  void* pointer2free);



/**
 * @fn void mem_backfreeall(int* memory_len_ptr,int* amount_elements_ptr,
 * 				  void*** all_mem_blocks_ptr,mem_exitcode_t* exitcode_ptr,
 * 				  int* showinfo_ptr);
 * @brief Frees all memory blocks (should never be called directly)
 *
 * Frees all memory blocks 
 * and deletes the array at "all_mem_blocks".
 *
 * @param memory_len_ptr Pointer to current len of "all_mem_blocks"\
							in the memory in units of sizeof(void*)
 * @param amount_elements_ptr Amount of memory blocks currently allocated
 * @param all_mem_blocks_ptr Array of all memory block pointers
 * @param exitcode_ptr Exitcode address (needed in later error processing)
 */
void mem_backfreeall(int* memory_len_ptr,int* amount_elements_ptr,
				  void*** all_mem_blocks_ptr,mem_exitcode_t* exitcode_ptr,
				  int* showinfo_ptr);



/**
 * @fn int mem_searcher(int* memory_len_ptr, void*** all_mem_blocks_ptr,\
						int* exitcode_ptr, void* pointer2search);
 * @brief Searches for the given pointer in all_mem_blocks\
				(should never be called directly)
 *
 * @param memory_len_ptr Pointer to current len of "all_mem_blocks"\
							in the memory in units of sizeof(void*)
 * @param all_mem_blocks_ptr Array of all memory block pointers
 * @param exitcode_ptr Exitcode address (needed in later error processing)
 * @param pointer2search Pointer to be searched for
 *
 * @return Last index in "all_mem_blocks", where "pointer2search" was found
 */
int mem_searcher(int* memory_len_ptr, void*** all_mem_blocks_ptr,
				mem_exitcode_t* exitcode_ptr, void* pointer2search);

#endif //ifndef GENERALMEMORYMANAGER_H
