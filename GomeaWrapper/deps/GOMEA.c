/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Header -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * GOMEA.c
 *
 * Copyright (c) Peter A.N. Bosman
 *
 * The software in this file is the proprietary information of
 * Peter A.N. Bosman.
 *
 * IN NO EVENT WILL THE AUTHOR OF THIS SOFTWARE BE LIABLE TO YOU FOR ANY
 * DAMAGES, INCLUDING BUT NOT LIMITED TO LOST PROFITS, LOST SAVINGS, OR OTHER
 * INCIDENTIAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR THE INABILITY
 * TO USE SUCH PROGRAM, EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGES, OR FOR ANY CLAIM BY ANY OTHER PARTY. THE AUTHOR MAKES NO
 * REPRESENTATIONS OR WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE
 * AUTHOR SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY ANYONE AS A RESULT OF
 * USING, MODIFYING OR DISTRIBUTING THIS SOFTWARE OR ITS DERIVATIVES.
 */
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=- Section Global Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
#if defined(_WIN32) || defined(WIN32) || defined(__CYGWIN__) || defined(__MINGW32__) || defined(__BORLANDC__)
#define OS_WIN
#endif
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include <stdio.h>
#include <stdlib.h>
#ifdef OS_WIN
#include <stdint.h>
#endif
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=- Section Global Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/

/* Special typedefs */
typedef void (*JuliaFunctionRunner)(void*, int, double*, double*, double*);
/* Global variables */
char       *terminated;                                             /* Whether a specific GOMEA with the restart scheme has terminated. */
double     *elitist_solution,                                       /* The very best solution ever evaluated. */
          **sostr,                                                  /* Set of solutions to reach. */
         ***populations,                                            /* The populations containing the solutions. */
         ***offsprings;                                             /* Offspring solutions (one set per population). */
short       write_generational_statistics,                          /* Whether to compute and write statistics every generation (0 = no). */
            write_generational_solutions,                           /* Whether to write the population every generation (0 = no). */
            print_verbose_overview,                                 /* Whether to print a overview of settings (0 = no). */
            print_FOSs_contents,                                    /* Whether to print the contents of the FOS structure each generation (0 = no). */
            use_ilse,                                               /* Whether incremental linkage subset evaluation is to be used. */
            vosostr_hit_status,                                     /* Whether the vosostr hit has happened yet: a solution has been evaluated and a value >= VTR or a STR has been found (0 = no, 1 = yes, 2 = yes, but this is no longer the first time. */
            vtr_exists,                                             /* Whether a vtr exists. */
            sostr_exists;                                           /* Whether a sostr exists. */
int         problem_index,                                          /* The index of the optimization problem. */
            FOSs_structure_index,                                   /* The index of the FOS structure. */
            number_of_parameters,                                   /* The number of parameters to be optimized. */
            number_of_solutions_in_sostr,                           /* The number of solutions in the set of solutions to reach. */
            number_of_generations,                                  /* The current generation count. */
           *population_sizes,                                       /* The number of solutions in each population. */
            base_population_size,                                   /* The minimum population size used in the smallest GOMEA instance. */
           *number_of_subgenerations_per_GOMEA,                     /* The number of subgenerations per GOMEA instance. */
            number_of_subgenerations_per_GOMEA_factor,              /* The factor by which the number of subgenerations increases with every new population. */
           *no_improvement_stretchs,                                /* The number of subsequent generations without an improvement for every GOMEA. */
            number_of_GOMEAs,                                       /* The number of GOMEAs currently running in multipop configuration. */
            maximum_number_of_GOMEAs,                               /* The maximum number of GOMEAs running in multipop configuration. */
         ***FOSs,                                                   /* The family of subsets linkage struture. */
          **FOSs_number_of_indices,                                 /* The number of variables in each linkage subset. */
           *FOSs_length,                                            /* The number of linkage subsets. */
            minimum_GOMEA_index;                                    /* The minimum GOMEA index that corresponds to the GOMEA that is still allowed to run (lower ones should be stopped because of average fitness being lower than that of a higher one). */
long        maximum_number_of_evaluations,                          /* The maximum number of evaluations. */
            maximum_number_of_milliseconds,                         /* The maximum number of milliseconds. */
            timestamp_start,                                        /* The time stamp in milliseconds for when the program was started. */
            timestamp_start_after_init,                             /* The time stamp in milliseconds for when the algorithm was started (after problem initialization). */
            number_of_evaluations,                                  /* The current number of times a function evaluation was performed. */
            elitist_solution_number_of_evaluations,                 /* The number of evaluations until the elitist solution. */
            elitist_solution_hitting_time,                          /* The hitting time of the elitist solution. */
            vosostr_number_of_evaluations,                          /* The number of evaluations until a solution that was to be reached (either vtr or in sostr). */
            vosostr_hitting_time;                                   /* The hitting time of a solution that was to be reached (either vtr or in sostr). */
long long   number_of_bit_flip_evaluations,                         /* The number of bit-flip evaluations. */
            elitist_solution_number_of_bit_flip_evaluations,        /* The number of bit-flip evaluations until the elitist solution. */
            vosostr_number_of_bit_flip_evaluations;                 /* The number of bit-flip evaluations until a solution that was to be reached (either vtr or in sostr). */
double      elitist_solution_objective_value,                       /* The objective value of the elitist solution. */
            elitist_solution_constraint_value,                      /* The constraint value of the elitist solution. */
            vtr,                                                    /* The value to reach (fitness of best solution that is feasible). */
          **objective_values,                                       /* Objective values for population members. */
          **constraint_values,                                      /* Sum of all constraint violations for population members. */
          **objective_values_offsprings,                            /* Objective values of selected solutions. */
          **constraint_values_offsprings,                           /* Sum of all constraint violations of selected solutions. */
           *objective_values_best_of_generation,
           *constraint_values_best_of_generation,
           *average_objective_values,
           *average_constraint_values,
         ***dependency_matrices;                                    /* Measure of dependency between any two variables (higher is more dependent). */
int64_t     random_seed,                                            /* The seed used for the random-number generator. */
            random_seed_changing;                                   /* Internally used variable for randomly setting a random seed. */

void *juliafunction;               /* Used to store the Julia function that should be evaluated. */
JuliaFunctionRunner calljulia; /* The function that calls the Julia function from a julia environment, 
                                      allowing for the use of julia like constructs such as closures */
char earlystop;
char reencode;
int randomkeysuntil;
char rescale;
char translate;
char eval_performs_ls;
char dependency_matrix_method;

/* Generate new solutions according to
  0 - Global [range_ub[0], range_ub[0]] range.
  1 - Per-index i a [rangeslb[i], rangesub[i]] range. */
char new_solution_generator;
double *rangeslb,
       *rangesub;
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
// TODO: Prune unused functions.
void *Malloc( long size );
int *mergeSortIntegersDecreasing( int *array, int array_size );
void mergeSortIntegersDecreasingWithinBounds( int *array, int *sorted, int *tosort, int p, int q );
void mergeSortIntegersDecreasingMerge( int *array, int *sorted, int *tosort, int p, int r, int q );
int *mergeSortDoublesIncreasing( double *array, int array_size );
void mergeSortDoublesIncreasingWithinBounds( double *array, int *sorted, int *tosort, int p, int q );
void mergeSortDoublesIncreasingMerge( double *array, int *sorted, int *tosort, int p, int r, int q );
void interpretCommandLine( int argc, char **argv );
void parseCommandLine( int argc, char **argv );
void parseOptions( int argc, char **argv, int *index );
void printAllInstalledProblems();
void printAllInstalledFOSs();
void optionError( char **argv, int index );
void parseParameters( int argc, char **argv, int *index );
void printUsage();
void checkOptions();
void printVerboseOverview();
double randomRealUniform01();
int randomInt( int maximum );
int *randomPermutation( int n );
char *installedProblemName( int index );
int numberOfInstalledProblems();
void installedProblemEvaluation( int index, double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *parameters_before, double objective_value_before, double constraint_value_before );
int *randomKeysToIntegerPermutation( double *parameters );
void sortFunctionProblemEvaluation( int *parameters_as_integer_permutation, double *objective_value, double *constraint_value );
void initializeNewGOMEA();
void initializeNewGOMEAMemory();
void initializeNewGOMEAPopulationAndFitnessValues();
void initializeValueAndSetOfSolutionsToReach();
short initializeValueToReach();
short initializeSetOfSolutionsToReach();
void initializeRandomNumberGenerator();
void selectFinalSurvivorsSpecificGOMEA( int gomea_index );
char betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y );
char equalFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y );
void writeGenerationalStatistics();
void writeGenerationalSolutions( char is_final_generation );
void writeRunningTime( char *filename );
void writeElitistSolution();
char checkTermination();
char checkNumberOfEvaluationsTerminationCondition();
char checkVOSOSTRTerminationCondition();
char checkNumberOfMilliSecondsTerminationCondition();
void generationalStepAllGOMEAs();
void makeOffspringSpecificGOMEA( int gomea_index );
char *installedFOSStructureName( int index );
int numberOfInstalledFOSStructures();
void learnFOSSpecificGOMEA( int gomea_index );
void learnUnivariateFOSSpecificGOMEA( int gomea_index );
int **learnLTFOSSpecificGOMEA( int gomea_index, short compute_dependency_matrices, short compute_parent_child_relations, int *number_of_parent_child_relations );
int determineNearestNeighbour( int index, double **S_matrix, int *mpm_number_of_indices, int mpm_length );
int **learnMLNFOSSpecificGOMEA( int gomea_index, short compute_dependency_matrices, short compute_parent_child_relations, int *number_of_parent_child_relations );
void learnLTNFOSSpecificGOMEA( int gomea_index );
void learnLTNFOSWithOrWithoutFilteringSpecificGOMEA( int gomea_index, short use_filtering );
void learnFilteredLTFOSSpecificGOMEA( int gomea_index, short compute_dependency_matrices );
void learnFilteredMLNFOSSpecificGOMEA( int gomea_index, short compute_dependency_matrices );
void learnFilteredLTNFOSSpecificGOMEA( int gomea_index );
void filterParentChildRelationsAndCreateNewFOSSpecificGOMEA( int gomea_index, int **parent_child_relations, int number_of_parent_child_relations );
double computeLinkageStrengthSpecificGOMEA( int gomea_index, int *variables, int number_of_variables );
void computeDependencyMatrixSpecificGOMEA( int gomea_index );
double *estimateParametersForSingleBinaryMarginalSpecificGOMEA( int gomea_index, int *indices, int number_of_indices, int *factor_size );
void uniquifyFOSSpecificGOMEA( int gomea_index );
short linkageSubsetsOfSameLengthAreDuplicates( int *linkageSubset0, int *linkageSubset1, int length );
void printFOSContentsSpecificGOMEA( int gomea_index );
double log2( double x );
void generateAndEvaluateNewSolutionsToFillOffspringSpecificGOMEA( int gomea_index );
double *generateAndEvaluateNewSolutionSpecificGOMEA( int gomea_index, int parent_index, double *obj, double *con );
void shuffleFOSSpecificGOMEA( int gomea_index );
void shuffleFOSSubsetsSpecificGOMEA( int gomea_index );
void ezilaitiniAllGOMEAs();
void ezilaitiniSpecificGOMEA( int gomea_index );
void ezilaitiniSpecificGOMEAMemoryForPopulationAndOffspring( int gomea_index );
void ezilaitiniValueAndSetOfSolutionsToReach();
void ezilaitiniProblem( int index );
long getMilliSecondsRunning();
long getMilliSecondsRunningAfterInit();
long getMilliSecondsRunningSinceTimeStamp( long timestamp );
long getCurrentTimeStampInMilliSeconds();
void run();
void multiPopGOMEA();
int main( int argc, char **argv );

#ifdef __cplusplus
extern "C" { 
#endif
// Entrypoint
// Follows most arguments from the GOMEA implementation 
// 1: problem (now a function, rather than integer)
// 2: number of parameters
// 3: fos index
// 4: maximum number of evaluations (-1 for no limit)
// 5: maximum number of milliseconds allowed (-1 for no limit)
void optimizeGOMEA(void *functrunk, JuliaFunctionRunner calljulia, 
  int dim, int fos, int eva, int mil, 
  int randomkeysuntila, 
  char reencodea, char rescalea, char translatea,
  char range_type, double *range_lb, double *range_ub,
  double *objective, double *constraint, double *best_solution,
  char eval_performs_ls_a,
  char dependency_matrix_method_a);
#ifdef __cplusplus
} 
#endif

#define MAX(a,b) \
  ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-= Section Elementary Operations -=-=-=-=-=-=-=-=-=-=-*/
/**
 * Allocates memory and exits the program in case of a memory allocation failure.
 */
void *Malloc( long size )
{
  void *result;

  result = (void *) malloc( size );
  if( !result )
  {
    printf( "\n" );
    printf( "Error while allocating memory in Malloc( %ld ), aborting program.", size );
    printf( "\n" );

    exit( 0 );
  }

  return( result );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Merge Sort -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Sorts an array of integers and returns the sort-order (large to small).
 */
int *mergeSortIntegersDecreasing( int *array, int array_size )
{
  int i, *sorted, *tosort;

  sorted = (int *) Malloc( array_size * sizeof( int ) );
  tosort = (int *) Malloc( array_size * sizeof( int ) );
  for( i = 0; i < array_size; i++ )
    tosort[i] = i;

  if( array_size == 1 )
    sorted[0] = 0;
  else
    mergeSortIntegersDecreasingWithinBounds( array, sorted, tosort, 0, array_size-1 );

  free( tosort );

  return( sorted );
}

/**
 * Subroutine of merge sort, sorts the part of the array between p and q.
 */
void mergeSortIntegersDecreasingWithinBounds( int *array, int *sorted, int *tosort, int p, int q )
{
  int r;

  if( p < q )
  {
    r = (p + q) / 2;
    mergeSortIntegersDecreasingWithinBounds( array, sorted, tosort, p, r );
    mergeSortIntegersDecreasingWithinBounds( array, sorted, tosort, r+1, q );
    mergeSortIntegersDecreasingMerge( array, sorted, tosort, p, r+1, q );
  }
}

/**
 * Subroutine of merge sort, merges the results of two sorted parts.
 */
void mergeSortIntegersDecreasingMerge( int *array, int *sorted, int *tosort, int p, int r, int q )
{
  int i, j, k, first;

  i = p;
  j = r;
  for( k = p; k <= q; k++ )
  {
    first = 0;
    if( j <= q )
    {
      if( i < r )
      {
        if( array[tosort[i]] > array[tosort[j]] )
          first = 1;
      }
    }
    else
      first = 1;

    if( first )
    {
      sorted[k] = tosort[i];
      i++;
    }
    else
    {
      sorted[k] = tosort[j];
      j++;
    }
  }

  for( k = p; k <= q; k++ )
    tosort[k] = sorted[k];
}

/**
 * Sorts an array of doubles and returns the sort-order (large to small).
 */
int *mergeSortDoublesIncreasing( double *array, int array_size )
{
  int i, *sorted, *tosort;

  sorted = (int *) Malloc( array_size * sizeof( int ) );
  tosort = (int *) Malloc( array_size * sizeof( int ) );
  for( i = 0; i < array_size; i++ )
    tosort[i] = i;

  if( array_size == 1 )
    sorted[0] = 0;
  else
    mergeSortDoublesIncreasingWithinBounds( array, sorted, tosort, 0, array_size-1 );

  free( tosort );

  return( sorted );
}

/**
 * Subroutine of merge sort, sorts the part of the array between p and q.
 */
void mergeSortDoublesIncreasingWithinBounds( double *array, int *sorted, int *tosort, int p, int q )
{
  int r;

  if( p < q )
  {
    r = (p + q) / 2;
    mergeSortDoublesIncreasingWithinBounds( array, sorted, tosort, p, r );
    mergeSortDoublesIncreasingWithinBounds( array, sorted, tosort, r+1, q );
    mergeSortDoublesIncreasingMerge( array, sorted, tosort, p, r+1, q );
  }
}

/**
 * Subroutine of merge sort, merges the results of two sorted parts.
 */
void mergeSortDoublesIncreasingMerge( double *array, int *sorted, int *tosort, int p, int r, int q )
{
  int i, j, k, first;

  i = p;
  j = r;
  for( k = p; k <= q; k++ )
  {
    first = 0;
    if( j <= q )
    {
      if( i < r )
      {
        if( array[tosort[i]] < array[tosort[j]] )
          first = 1;
      }
    }
    else
      first = 1;

    if( first )
    {
      sorted[k] = tosort[i];
      i++;
    }
    else
    {
      sorted[k] = tosort[j];
      j++;
    }
  }

  for( k = p; k <= q; k++ )
    tosort[k] = sorted[k];
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/




/*-=-=-=-=-=-=-=-=-=-=- Section Interpret Command Line -=-=-=-=-=-=-=-=-=-=-*/
/**
 * Parses and checks the command line.
 */
void interpretCommandLine( int argc, char **argv )
{
  parseCommandLine( argc, argv );
  
  checkOptions();
}

/**
 * Parses the command line.
 * For options, see printUsage.
 */
void parseCommandLine( int argc, char **argv )
{
  int index;

  index = 1;

  parseOptions( argc, argv, &index );
  
  parseParameters( argc, argv, &index );
}

/**
 * Parses only the options from the command line.
 */
void parseOptions( int argc, char **argv, int *index )
{
  double dummy;

  write_generational_statistics = 0;
  write_generational_solutions  = 0;
  print_verbose_overview        = 0;
  print_FOSs_contents           = 0;
  use_ilse                      = 0;

  for( ; (*index) < argc; (*index)++ )
  {
    if( argv[*index][0] == '-' )
    {
      /* If it is a negative number, the option part is over */
      if( sscanf( argv[*index], "%lf", &dummy ) && argv[*index][1] != '\0' )
        break;

      if( argv[*index][1] == '\0' )
        optionError( argv, *index );
      else if( argv[*index][2] != '\0' )
        optionError( argv, *index );
      else
      {
        switch( argv[*index][1] )
        {
          case '?': printUsage(); break;
          case 'P': printAllInstalledProblems(); break;
          case 'F': printAllInstalledFOSs(); break;
          case 's': write_generational_statistics = 1; break;
          case 'w': write_generational_solutions  = 1; break;
          case 'v': print_verbose_overview        = 1; break;
          case 'f': print_FOSs_contents           = 1; break;
          case 'i': use_ilse                      = 1; break;
          default : optionError( argv, *index );
        }
      }
    }
    else /* Argument is not an option, so option part is over */
     break;
  }
}

/**
 * Writes the names of all installed problems to the standard output.
 */
void printAllInstalledProblems()
{
  int i, n;

  n = numberOfInstalledProblems();
  printf( "Installed optimization problems:\n" );
  for( i = 0; i < n; i++ )
    printf( "%3d: %s\n", i, installedProblemName( i ) );

  // exit( 0 );
}

/**
 * Writes the names of all installed FOS structures to the standard output.
 */
void printAllInstalledFOSs()
{
  int i, n;

  n = numberOfInstalledFOSStructures();
  printf( "Installed FOS structures:\n" );
  for( i = 0; i < n; i++ )
    printf( "%3d: %s\n", i, installedFOSStructureName( i ) );

  // exit( 0 );
}

/**
 * Informs the user of an illegal option and exits the program.
 */
void optionError( char **argv, int index )
{
  printf( "Illegal option: %s\n\n", argv[index] );

  printUsage();
}

/**
 * Parses only the EA parameters from the command line.
 */
void parseParameters( int argc, char **argv, int *index )
{
  int noError;

  if( (argc - *index) != 5 )
  {
    printf( "Number of parameters is incorrect, require 5 parameters (you provided %d).\n\n", (argc - *index) );

    printUsage();
  }

  noError = 1;
  noError = noError && sscanf( argv[*index+0], "%d", &problem_index );
  noError = noError && sscanf( argv[*index+1], "%d", &number_of_parameters );
  noError = noError && sscanf( argv[*index+2], "%d", &FOSs_structure_index );
  noError = noError && sscanf( argv[*index+3], "%ld", &maximum_number_of_evaluations );
  noError = noError && sscanf( argv[*index+4], "%ld", &maximum_number_of_milliseconds );

  if( !noError )
  {
    printf( "Error parsing parameters.\n\n" );

    printUsage();
  }
}

/**
 * Prints usage information and exits the program.
 */
void printUsage()
{
  printf("Usage: GOMEA [-?] [-P] [-F] [-s] [-w] [-v] [-f] [-i] pro dim fos eva mil\n");
  printf("   -?: Prints out this usage information.\n");
  printf("   -P: Prints out a list of all installed optimization problems.\n");
  printf("   -F: Prints out a list of all installed FOS structures.\n");
  printf("   -s: Enables computing and writing of statistics every generation.\n");
  printf("   -w: Enables writing of solutions and their fitnesses every generation.\n");
  printf("   -v: Enables verbose mode. Prints the settings before starting the run.\n");
  printf("   -f: Enables printing the contents of the FOS every generation.\n");
  printf("   -i: Enables Incremental Linkage Subset Evaluation (ILSE).\n");
  printf("\n");
  printf("  pro: Index of optimization problem to be solved (minimization).\n");
  printf("  dim: Number of parameters.\n");
  printf("  fos: Index of FOS structure to use.\n");
  printf("  eva: Maximum number of evaluations allowed (-1 for no limit).\n");
  printf("  mil: Maximum number of milliseconds allowed (-1 for no limit).\n");
  // exit( 0 );
}

/**
 * Checks whether the selected options are feasible.
 */
void checkOptions()
{
  if( number_of_parameters < 1 )
  {
    printf("\n");
    printf("Error: number of parameters < 1 (read: %d). Require number of parameters >= 1.", number_of_parameters);
    printf("\n\n");

    // exit( 0 );
  }

  if( installedProblemName( problem_index ) == NULL )
  {
    printf("\n");
    printf("Error: unknown index for problem (read index %d).", problem_index );
    printf("\n\n");

    // exit( 0 );
  }

  if( installedFOSStructureName( FOSs_structure_index ) == NULL )
  {
    printf("\n");
    printf("Error: unknown index for FOS structure (read index %d).", FOSs_structure_index );
    printf("\n\n");

    // exit( 0 );
  }
}

/**
 * Prints the settings as read from the command line.
 */
void printVerboseOverview()
{
  printf("### Settings ######################################\n");
  printf("#\n");
  printf("# Statistics writing every generation: %s\n", write_generational_statistics ? "enabled" : "disabled");
  printf("# Population file writing            : %s\n", write_generational_solutions ? "enabled" : "disabled");
  printf("# FOS contents printing              : %s\n", print_FOSs_contents ? "enabled" : "disabled");
  printf("# ILSE                               : %s\n", use_ilse ? "enabled" : "disabled" );
  printf("#\n");
  printf("###################################################\n");
  printf("#\n");
  printf("# Problem                            = %s\n", installedProblemName( problem_index ));
  printf("# Number of parameters               = %d\n", number_of_parameters);
  printf("# FOS structure                      = %s\n", installedFOSStructureName( FOSs_structure_index ));
  printf("# Maximum number of evaluations      = %ld\n", maximum_number_of_evaluations);
  printf("# Maximum number of milliseconds     = %ld\n", maximum_number_of_milliseconds);
  printf("# Random seed                        = %ld\n", (long) random_seed);
  printf("#\n");
  printf("###################################################\n");
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Random Numbers -=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns a random double, distributed uniformly between 0 and 1.
 */
double randomRealUniform01()
{
  int64_t n26, n27;
  double  result;

  random_seed_changing = (random_seed_changing * 0x5DEECE66DLLU + 0xBLLU) & ((1LLU << 48) - 1);
  n26                  = (int64_t)(random_seed_changing >> (48 - 26));
  random_seed_changing = (random_seed_changing * 0x5DEECE66DLLU + 0xBLLU) & ((1LLU << 48) - 1);
  n27                  = (int64_t)(random_seed_changing >> (48 - 27));
  result               = (((int64_t)n26 << 27) + n27) / ((double) (1LLU << 53));

  return( result );
}

double randomRealUniformAB(double a, double b)
{
  double result;

  result = randomRealUniform01() * (b - a) + a;

  return( result );
}
        
/**
 * Returns a random integer, distributed uniformly between 0 and maximum.
 */
int randomInt( int maximum )
{
  int result;
  
  result = (int) (((double) maximum)*randomRealUniform01());
  
  return( result );
}

/**
 * Returns a random compact (using integers 0,1,...,n-1) permutation
 * of length n using the Fisher-Yates shuffle.
 */
int *randomPermutation( int n )
{
  int i, j, dummy, *result;

  result = (int *) Malloc( n*sizeof( int ) );
  for( i = 0; i < n; i++ )
    result[i] = i;

  for( i = n-1; i > 0; i-- )
  {
    j         = randomInt( i+1 );
    dummy     = result[j];
    result[j] = result[i];
    result[i] = dummy;
  }

  return( result );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Problems -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns the name of an installed problem.
 */
char *installedProblemName( int index )
{
  switch( index )
  {
    case  0: return( (char *) "Provided Function" );                 
  }

  return( NULL );
}

/**
 * Returns the number of problems installed.
 */
int numberOfInstalledProblems()
{
  static int result = -1;
  
  if( result == -1 )
  {
    result = 0;
    while( installedProblemName( result ) != NULL )
      result++;
  }
  
  return( result );
}

/**
 * Computes the value of the single objective
 * and the sum of all constraint violations
 * function.
 */
void installedProblemEvaluation( int index, double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *parameters_before, double objective_value_before, double constraint_value_before )
{
  short same;
  int   i, j;

  /* Count the evaluation */
  number_of_evaluations++;
  if( number_of_touched_parameters == 0 )
    number_of_bit_flip_evaluations += number_of_parameters;
  else
  {
    for( i = 0; i < number_of_touched_parameters; i++ )
    {
      if( parameters[touched_parameters_indices[i]] != parameters_before[touched_parameters_indices[i]] )
        number_of_bit_flip_evaluations++;
    }
  }

  /* Do the actual evaluation */
  *objective_value  = 0.0;
  *constraint_value = 0.0;
  
    // if ( earlystop )
    //   printf("Notice: evaluation after stopping criterion was hit.");

  // Call the function pointer to call the Julia function to optimize in the functrunk.
  (*calljulia)(juliafunction, number_of_parameters, parameters, objective_value, constraint_value);

  /* Check the VTR */
  if( !vosostr_hit_status )
  {
    if( vtr_exists > 0 )
    {
      if( ((*constraint_value) == 0) && ((*objective_value) >= vtr)  )
      {
        vosostr_hit_status = 1;
      }
    }
  }

  /* Check the SOSTR */
  if( !vosostr_hit_status )
  {
    if( sostr_exists )
    {
      for( i = 0; i < number_of_solutions_in_sostr; i++ )
      {
        same = 1;
        for( j = 0; j < number_of_parameters; j++ )
        {
          if( parameters[j] != sostr[i][j] )
          {
            same = 0;
            break;
          }
        }
        if( same )
        {
          vosostr_hit_status = 1;
          break;
        }
      }
    }
  }

  /* Check the VOSOSTR */
  if( vosostr_hit_status == 1 )
  {
    vosostr_hit_status                     = 2;
    vosostr_hitting_time                   = getMilliSecondsRunningAfterInit();
    vosostr_number_of_evaluations          = number_of_evaluations;
    vosostr_number_of_bit_flip_evaluations = number_of_bit_flip_evaluations;

    // writeRunningTime( (char *) "vosostr_hitting_time.dat" );
  }

  /* Update elitist solution */
  if( (number_of_evaluations == 1) || betterFitness( *objective_value, *constraint_value, elitist_solution_objective_value, elitist_solution_constraint_value ) )
  {
    for( i = 0; i < number_of_parameters; i++ )
      elitist_solution[i] = parameters[i];

    elitist_solution_objective_value                = *objective_value;
    elitist_solution_constraint_value               = *constraint_value;
    elitist_solution_hitting_time                   = getMilliSecondsRunningAfterInit();
    elitist_solution_number_of_evaluations          = number_of_evaluations;
    elitist_solution_number_of_bit_flip_evaluations = number_of_bit_flip_evaluations;

    // writeRunningTime( (char *) "elitist_solution_hitting_time.dat" );

    writeElitistSolution();
  }

  /* Exit early, depending on VOSOSTR status */
  if( vosostr_hit_status != 0 ) {
    earlystop = 1;
    // Not using the vosostr. So I doubt this has any impact.
  }
  if( checkNumberOfEvaluationsTerminationCondition() ) {
    earlystop = 1;
    // This is already handled in the termination criterion.
  }
  if( checkNumberOfMilliSecondsTerminationCondition() ) {
    earlystop = 1;
    // This is already handled in the termination criterion.
  }
}

int *randomKeysToIntegerPermutation( double *parameters )
{
  int *result;

  result = mergeSortDoublesIncreasing( parameters, number_of_parameters );
  
  return( result );
}

void sortFunctionProblemEvaluation( int *parameters, double *objective_value, double *constraint_value )
{
  int    i, j;
  double result;

  result = 0.0;
  for( i = 0; i < number_of_parameters; i++ )
    for( j = i+1; j < number_of_parameters; j++ )
      result += parameters[i] > parameters[j] ? 0 : 1;

  *objective_value  = result;
  *constraint_value = 0;
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Initialize -=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * Performs initialization for a single GOMEA.
 */
void initializeNewGOMEA()
{
  int gomea_index;

  if( number_of_GOMEAs == 0 )
  {
    populations                  = (double ***) Malloc( maximum_number_of_GOMEAs*sizeof( double ** ) );
    objective_values             = (double **) Malloc( maximum_number_of_GOMEAs*sizeof( double * ) );
    constraint_values            = (double **) Malloc( maximum_number_of_GOMEAs*sizeof( double * ) );
    offsprings                   = (double ***) Malloc( maximum_number_of_GOMEAs*sizeof( double ** ) );
    objective_values_offsprings  = (double **) Malloc( maximum_number_of_GOMEAs*sizeof( double * ) );
    constraint_values_offsprings = (double **) Malloc( maximum_number_of_GOMEAs*sizeof( double * ) );
    average_objective_values     = (double *) Malloc( maximum_number_of_GOMEAs*sizeof( double ) );
    average_constraint_values    = (double *) Malloc( maximum_number_of_GOMEAs*sizeof( double ) );
    terminated                   = (char *) Malloc( maximum_number_of_GOMEAs*sizeof( char ) );
    dependency_matrices          = (double ***) Malloc( maximum_number_of_GOMEAs*sizeof( double ** ) );
    FOSs                         = (int ***) Malloc( maximum_number_of_GOMEAs*sizeof( int ** ) );
    FOSs_number_of_indices       = (int **) Malloc( maximum_number_of_GOMEAs*sizeof( int * ) );
    FOSs_length                  = (int *) Malloc( maximum_number_of_GOMEAs*sizeof( int ) );
    objective_values_best_of_generation  = (double *) Malloc( maximum_number_of_GOMEAs*sizeof( double ) );
    constraint_values_best_of_generation = (double *) Malloc( maximum_number_of_GOMEAs*sizeof( double ) );

    population_sizes[number_of_GOMEAs]                   = base_population_size;
    number_of_subgenerations_per_GOMEA[number_of_GOMEAs] = 1;
  }
  else
  {
    population_sizes[number_of_GOMEAs]                   = 2*population_sizes[number_of_GOMEAs-1];
    number_of_subgenerations_per_GOMEA[number_of_GOMEAs] = 1;
    for( gomea_index = 0; gomea_index < number_of_GOMEAs; gomea_index++ )
      number_of_subgenerations_per_GOMEA[gomea_index] *= number_of_subgenerations_per_GOMEA_factor;
  }

  terminated[number_of_GOMEAs]              = 0;
  no_improvement_stretchs[number_of_GOMEAs] = 0;
  FOSs[number_of_GOMEAs]                    = NULL;

  initializeNewGOMEAMemory();

  initializeNewGOMEAPopulationAndFitnessValues();

  number_of_GOMEAs++;
}

/**
 * Initializes the memory for a single GOMEA.
 */
void initializeNewGOMEAMemory()
{
  int i;

  populations[number_of_GOMEAs]                  = (double **) Malloc( population_sizes[number_of_GOMEAs]*sizeof( double * ) );
  objective_values[number_of_GOMEAs]             = (double *) Malloc( population_sizes[number_of_GOMEAs]*sizeof( double ) );
  constraint_values[number_of_GOMEAs]            = (double *) Malloc( population_sizes[number_of_GOMEAs]*sizeof( double ) );
  offsprings[number_of_GOMEAs]                   = (double **) Malloc( population_sizes[number_of_GOMEAs]*sizeof( double * ) );
  objective_values_offsprings[number_of_GOMEAs]  = (double *) Malloc( population_sizes[number_of_GOMEAs]*sizeof( double ) );
  constraint_values_offsprings[number_of_GOMEAs] = (double *) Malloc( population_sizes[number_of_GOMEAs]*sizeof( double ) );

  for( i = 0; i < population_sizes[number_of_GOMEAs]; i++ )
    populations[number_of_GOMEAs][i] = (double *) Malloc( number_of_parameters*sizeof( double ) );

  for( i = 0; i < population_sizes[number_of_GOMEAs]; i++ )
    offsprings[number_of_GOMEAs][i] = (double *) Malloc( number_of_parameters*sizeof( double ) );

  dependency_matrices[number_of_GOMEAs] = (double **) Malloc( number_of_parameters*sizeof( double * ) );
  for( i = 0; i < number_of_parameters; i++ )
    (dependency_matrices[number_of_GOMEAs])[i] = (double *) Malloc( number_of_parameters*sizeof( double ) );

  FOSs[number_of_GOMEAs] = NULL;
}

/**
 * Initializes the population and the objective values by randomly
 * generation n solutions.
 */
void initializeNewGOMEAPopulationAndFitnessValues()
{
  int    i, j;
  double obj, con;

  objective_values_best_of_generation[number_of_GOMEAs]  = -1e+308;
  constraint_values_best_of_generation[number_of_GOMEAs] = 1e+308;
  for( i = 0; i < population_sizes[number_of_GOMEAs]; i++ )
  {
    // EDIT: Here a new random solution is generated.
    // Value now ranges from rangeslb[j * new_solution_generator]
    // Until rangesub[j * new_solution_generator]
    for( j = 0; j < number_of_parameters; j++ )
      if( new_solution_generator == 0 )
        populations[number_of_GOMEAs][i][j] = randomRealUniformAB(rangeslb[0], rangesub[0]);
      else
        populations[number_of_GOMEAs][i][j] = randomRealUniformAB(rangeslb[j], rangesub[j]);

    installedProblemEvaluation( problem_index, populations[number_of_GOMEAs][i], &obj, &con, 0, NULL, NULL, 0, 0 );
    objective_values[number_of_GOMEAs][i]  = obj;
    constraint_values[number_of_GOMEAs][i] = con;

    if( betterFitness( objective_values[number_of_GOMEAs][i], constraint_values[number_of_GOMEAs][i], objective_values_best_of_generation[number_of_GOMEAs], constraint_values_best_of_generation[number_of_GOMEAs] ) )
    {
      objective_values_best_of_generation[number_of_GOMEAs]  = objective_values[number_of_GOMEAs][i];
      constraint_values_best_of_generation[number_of_GOMEAs] = constraint_values[number_of_GOMEAs][i];
    }
  }
}

/**
 * Checks to see if files exists with values and solutions to reach.
 */
void initializeValueAndSetOfSolutionsToReach()
{
  vtr_exists = 0;
  if( initializeValueToReach() )
    vtr_exists = 1;

  sostr_exists = 0;
  if( initializeSetOfSolutionsToReach() )
    sostr_exists = 1;
}

/**
 * Attempts to read the value to reach.
 */
short initializeValueToReach()
{
  char  c, string[100000], filename[1000];
  int   i;
  FILE *file;

  sprintf( filename, "vtr.txt" );
  file = fopen( filename, "r" );
  if( file == NULL )
    return( 0 );

  i    = 0;
  c    = fgetc( file );
  while( (c != '\n') && (c != EOF) && (c != ' ') )
  {
    string[i] = (char) c;
    c         = fgetc( file );
    i++;
  }
  string[i] = '\0';
  sscanf( string, "%le", &vtr);

  fclose( file );

  return( 1 );
}

/**
 * Attempts to read assets of solutions to reach.
 */
short initializeSetOfSolutionsToReach()
{
  char  c, string[100000], filename[1000];
  int   i, j;
  FILE *file;

  number_of_solutions_in_sostr = 0;

  sprintf( filename, "sostr.txt" );
  file = fopen( filename, "r" );
  if( file == NULL )
    return( 0 );

  do
  {
    c = fgetc( file );
    if( (c == '0') || (c == '1') )
    {
      number_of_solutions_in_sostr++;
      do
      {
        c = fgetc( file );
      } while( (c == '0') || (c == '1') );
    }
  }
  while( c != EOF );
  fclose( file );

  sostr = (double **) Malloc( number_of_solutions_in_sostr*sizeof( double * ) );
  for( i = 0; i < number_of_solutions_in_sostr; i++ )
    sostr[i] = (double *) Malloc( number_of_parameters*sizeof( double ) );

  file = fopen( filename, "r" );
  i    = 0;
  do
  {
    c = fgetc( file );
    if( c == EOF )
      break;
    if( !((c == '0') || (c == '1')) )
      continue;
    j = 0;
    while( (c != '\n') && (c != EOF) && (c != ' ') )
    {
      string[j] = (char) c;
      c         = fgetc( file );
      j++;
    }
    if( j != number_of_parameters )
    {
      printf("Error while reading %s: the number of parameters in at least one of the solutions does not match the number of parameters on the commandline (%d != %d)\n",filename,j,number_of_parameters);
      exit(0);
    }
    for( j = 0; j < number_of_parameters; j++ )
      sostr[i][j] = string[j] == '0' ? 0 : 1;
  }
  while( c != EOF );

  fclose( file );

  return( 1 );
}

/**
 * Initializes the pseudo-random number generator.
 */
void initializeRandomNumberGenerator()
{
  struct timeval tv;
  struct tm *timep;
  time_t tvwin;

  while( random_seed_changing == 0 )
  {
    #ifdef OS_WIN
    time(&tvwin);
    timep = localtime (&tvwin);
    #else
    gettimeofday( &tv, NULL );
    timep = localtime (&tv.tv_sec);
    #endif
    random_seed_changing = timep->tm_hour * 3600 * 1000 + timep->tm_min * 60 * 1000 + timep->tm_sec * 1000 + tv.tv_usec / 1000;
  }

  random_seed = random_seed_changing;
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*=-=-=-=-=-=-=-=-=-=-= Section Survivor Selection =-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Determines the solutions that finally survive the generation (offspring only).
 */
void selectFinalSurvivorsSpecificGOMEA( int gomea_index )
{
  int    i, j;
  double objective_values_best_of_generation_before, constraint_values_best_of_generation_before;

  objective_values_best_of_generation_before  = objective_values_best_of_generation[gomea_index];
  constraint_values_best_of_generation_before = constraint_values_best_of_generation[gomea_index];

  for( i = 0; i < population_sizes[gomea_index]; i++ )
  {
    for( j = 0; j < number_of_parameters; j++ )
      populations[gomea_index][i][j] = offsprings[gomea_index][i][j];
    objective_values[gomea_index][i]  = objective_values_offsprings[gomea_index][i];
    constraint_values[gomea_index][i] = constraint_values_offsprings[gomea_index][i];

    if( betterFitness( objective_values[gomea_index][i], constraint_values[gomea_index][i], objective_values_best_of_generation[gomea_index], constraint_values_best_of_generation[gomea_index] ) )
    {
      objective_values_best_of_generation[gomea_index]  = objective_values[gomea_index][i];
      constraint_values_best_of_generation[gomea_index] = constraint_values[gomea_index][i];
    }
  }

  if( !betterFitness( objective_values_best_of_generation[gomea_index], constraint_values_best_of_generation[gomea_index], objective_values_best_of_generation_before, constraint_values_best_of_generation_before ) )
    no_improvement_stretchs[gomea_index]++;
  else
    no_improvement_stretchs[gomea_index] = 0;
}

/**
 * Returns 1 if x is better than y, 0 otherwise.
 * x is not better than y unless:
 * - x and y are both infeasible and x has a smaller sum of constraint violations, or
 * - x is feasible and y is not, or
 * - x and y are both feasible and x has a larger objective value than y
 */
char betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y )
{
  char result;
  
  result = 0;

  if( constraint_value_x > 0 ) /* x is infeasible */
  {
    if( constraint_value_y > 0 ) /* Both are infeasible */
    {
      if( constraint_value_x < constraint_value_y )
       result = 1;
    }
  }
  else /* x is feasible */
  {
    if( constraint_value_y > 0 ) /* x is feasible and y is not */
      result = 1;
    else /* Both are feasible */
    {
      if( objective_value_x > objective_value_y )
        result = 1;
    }
  }

  return( result );
}

/**
 * Returns 1 if x is equally preferable to y, 0 otherwise.
 */
char equalFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y )
{
  char result;
  
  result = 0;

  if( (constraint_value_x == constraint_value_y) && (objective_value_x == objective_value_y) )
    result = 1;

  return( result );
}

void computeAverageFitnessSpecificGOMEA( int GOMEA_index )
{
  int i;

  average_objective_values[GOMEA_index]  = 0;
  average_constraint_values[GOMEA_index] = 0;
  for( i = 0; i < population_sizes[GOMEA_index]; i++ )
  {
    average_objective_values[GOMEA_index]  += objective_values[GOMEA_index][i];
    average_constraint_values[GOMEA_index] += constraint_values[GOMEA_index][i];
  }
  average_objective_values[GOMEA_index] /= (double) (population_sizes[GOMEA_index]);
  average_constraint_values[GOMEA_index] /= (double) (population_sizes[GOMEA_index]);
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Output =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Writes (appends) statistics about the current generation to a
 * file named "statistics.dat".
 */
void writeGenerationalStatistics()
{
  char    filename[1000], string[10000];
  int     i, n, gomea_index;
  double  objective_avg, objective_var, objective_best, objective_worst,
          constraint_avg, constraint_var, constraint_best, constraint_worst;
  FILE   *file;

  /* First compute the statistics */
  /* Average, best and worst */ 
  objective_avg    = 0.0;
  constraint_avg   = 0.0;
  objective_best   = objective_values[0][0];
  objective_worst  = objective_values[0][0];
  constraint_best  = constraint_values[0][0];
  constraint_worst = constraint_values[0][0];
  n                = 0;
  for( gomea_index = 0; gomea_index < number_of_GOMEAs; gomea_index++ )
  {
    for( i = 0; i < population_sizes[gomea_index]; i++ )
    {
      objective_avg += objective_values[gomea_index][i];
      constraint_avg += constraint_values[gomea_index][i];
      if( betterFitness( objective_values[gomea_index][i], constraint_values[gomea_index][i], objective_best, constraint_best ) )
      {
        objective_best  = objective_values[gomea_index][i];
        constraint_best = constraint_values[gomea_index][i];
      }
      if( betterFitness( objective_worst, constraint_worst, objective_values[gomea_index][i], constraint_values[gomea_index][i] ) )
      {
        objective_worst  = objective_values[gomea_index][i];
        constraint_worst = constraint_values[gomea_index][i];
      }
      n++;
    }
  }
  objective_avg = objective_avg / ((double) n);
  constraint_avg = constraint_avg / ((double) n);

  /* Variance */
  objective_var  = 0.0;
  constraint_var = 0.0;
  for( gomea_index = 0; gomea_index < number_of_GOMEAs; gomea_index++ )
  {
    for( i = 0; i < population_sizes[gomea_index]; i++ )
    {
      objective_var += (objective_values[gomea_index][i] - objective_avg)*(objective_values[gomea_index][i] - objective_avg);
      constraint_var += (constraint_values[gomea_index][i] - constraint_avg)*(constraint_values[gomea_index][i] - constraint_avg);
    }
  }
  objective_var = objective_var / ((double) n);
  constraint_var = constraint_var / ((double) n);

  if( objective_var <= 0.0 )
     objective_var = 0.0;
  if( constraint_var <= 0.0 )
     constraint_var = 0.0;

  /* Then write them */
  sprintf( filename, "statistics.dat" );
  file = NULL;
  if( number_of_generations == 0 )
  {
    file = fopen( filename, "w" );

    sprintf( string, "# Generation Evaluations      Average-obj.     Variance-obj.         Best-obj.        Worst-obj.        Elite-obj.      Average-con.     Variance-con.         Best-con.        Worst-con.        Elite-con.\n");
    fputs( string, file );
  }
  else
    file = fopen( filename, "a" );

  sprintf( string, "  %10d %11ld %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e\n", number_of_generations, number_of_evaluations, objective_avg, objective_var, objective_best, objective_worst, elitist_solution_objective_value, constraint_avg, constraint_var, constraint_best, constraint_worst, elitist_solution_constraint_value );
  fputs( string, file );

  fclose( file );
}

/**
 * Writes the solutions to various files. The filenames
 * contain the generation. If the flag is_final_generation
 * is set the generation number in the filename
 * is replaced with the word "final".
 *
 * populations_xxxxx_generation_xxxxx.dat: the populations
 * offsprings_xxxxx_generation_xxxxx.dat : the offsprings
 */
void writeGenerationalSolutions( char is_final_generation )
{
  char  filename[1000], string[10000];
  int   i, j, gomea_index;
  FILE *file;

  /* Populations */
  if( is_final_generation )
    sprintf( filename, "populations_generation_final.dat" );
  else
    sprintf( filename, "populations_generation_%05d.dat", number_of_generations );
  file = fopen( filename, "w" );

  for( gomea_index = 0; gomea_index < number_of_GOMEAs; gomea_index++ )
  {
    for( i = 0; i < population_sizes[gomea_index]; i++ )
    {
      for( j = 0; j < number_of_parameters; j++ )
      {
        sprintf( string, "%lf ", populations[gomea_index][i][j] );
        fputs( string, file );
      }
      sprintf( string, "     %17.10e %17.10e\n", objective_values[gomea_index][i], constraint_values[gomea_index][i] );
      fputs( string, file );
    }
  }
  fclose( file );
  
  /* Offsprings */
  if( (number_of_generations > 0) && (!is_final_generation) )
  {
    sprintf( filename, "offsprings_generation_%05d.dat", number_of_generations-1 );
    file = fopen( filename, "w" );

    for( gomea_index = 0; gomea_index < number_of_GOMEAs; gomea_index++ )
    {
      for( i = 0; i < population_sizes[gomea_index]; i++ )
      {
        for( j = 0; j < number_of_parameters; j++ )
        {
          sprintf( string, "%lf ", offsprings[gomea_index][i][j] );
          fputs( string, file );
        }
        sprintf( string, "     %17.10e %17.10e\n", objective_values_offsprings[gomea_index][i], constraint_values_offsprings[gomea_index][i] );
        fputs( string, file );
      }
    }

    fclose( file );
  }
}

void writeRunningTime( char *filename )
{
  char  string[10000];
  FILE *file;

  file = fopen( filename, "w" );
  sprintf( string, "# Column 1: Total number of milliseconds.\n");
  fputs( string, file );
  sprintf( string, "# Column 2: Total number of milliseconds after initialization.\n");
  fputs( string, file );
  sprintf( string, "# Column 3: Total number of evaluations.\n"); 
  fputs( string, file );
  sprintf( string, "# Column 4: Total number of bit-flip evaluations.\n"); 
  fputs( string, file );
#ifdef OS_WIN
  sprintf( string, "%ld %ld %ld %I64d\n", getMilliSecondsRunning(), getMilliSecondsRunningAfterInit(), number_of_evaluations, number_of_bit_flip_evaluations );
#else
  sprintf( string, "%ld %ld %ld %lld\n", getMilliSecondsRunning(), getMilliSecondsRunningAfterInit(), number_of_evaluations, number_of_bit_flip_evaluations );
#endif
  fputs( string, file );
  fclose( file );
}

void writeElitistSolution()
{
  
  // char       string[10000];
  // int        i, *elitist_solution_as_integer_permutation;
  // FILE      *file;

  // elitist_solution_as_integer_permutation = randomKeysToIntegerPermutation( elitist_solution );
  // file = fopen( "elitist_solution.dat", "w" );
  // sprintf( string, "# Column 1: Solution.\n");
  // fputs( string, file );
  // sprintf( string, "# Column 2: Objective value.\n");
  // fputs( string, file );
  // sprintf( string, "# Column 3: Constraint value.\n");
  // fputs( string, file );
  // for( i = 0; i < number_of_parameters; i++ )
  // {
  //   sprintf( string, "%d ", elitist_solution_as_integer_permutation[i] );
  //   fputs( string, file );
  // }
  // sprintf( string, "     %17.10e %17.10e\n", elitist_solution_objective_value, elitist_solution_constraint_value );
  // fputs( string, file );
  // fclose( file );
  // free( elitist_solution_as_integer_permutation );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Termination -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns 1 if termination should be enforced at the end of a generation, 0 otherwise.
 */
char checkTermination()
{
  if( checkNumberOfEvaluationsTerminationCondition() )
    return( 1 );
  
  if( checkNumberOfMilliSecondsTerminationCondition() )
    return( 1 );

  if( checkVOSOSTRTerminationCondition() )
      return( 1 );

  return( 0 );
}

/**
 * Returns 1 if the maximum number of evaluations
 * has been reached, 0 otherwise.
 */
char checkNumberOfEvaluationsTerminationCondition()
{
  if( (maximum_number_of_evaluations >= 0) && (number_of_evaluations >= maximum_number_of_evaluations) )
    return( 1 );

  return( 0 );
}

/**
 * Returns 1 if the value-to-reach has been reached.
 */
char checkVOSOSTRTerminationCondition()
{
  if( vosostr_hit_status > 0 )
    return( 1 );

  return( 0 );
}

/**
 * Returns 1 if the maximum number of milliseconds
 * has passed, 0 otherwise.
 */
char checkNumberOfMilliSecondsTerminationCondition()
{
  if( (maximum_number_of_milliseconds >= 0) && (getMilliSecondsRunning() > maximum_number_of_milliseconds) )
    return( 1 );

  return( 0 );
}

char checkTerminationSpecificGOMEA( int GOMEA_index )
{
  int i, j;

  // EDIT: NOT PREMATURELY STOPPING SMALLER POPULATIONS
/*
  for( i = GOMEA_index+1; i < number_of_GOMEAs; i++ )
  {
    if( betterFitness( average_objective_values[i], average_constraint_values[i], average_objective_values[GOMEA_index], average_constraint_values[GOMEA_index] ) )
    {
      minimum_GOMEA_index = GOMEA_index+1;

      return( 1 );
    }
  }
*/

  for( i = 1; i < population_sizes[GOMEA_index]; i++ )
  {
    for( j = 0; j < number_of_parameters; j++ )
    {
      if( populations[GOMEA_index][i][j] != populations[GOMEA_index][0][j] )
        return( 0 );
    }
  }

  return( 1 );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Variation -==-=-=-=-=-=-=-=-=-=-=-=-=-=*/
void generationalStepAllGOMEAsRecursiveFold( int GOMEA_index_smallest, int GOMEA_index_biggest );
void generationalStepAllGOMEAs()
{
  int GOMEA_index_smallest, GOMEA_index_biggest;

  GOMEA_index_biggest  = number_of_GOMEAs-1;
  GOMEA_index_smallest = 0;
  while( GOMEA_index_smallest <= GOMEA_index_biggest )
  {
    if( !terminated[GOMEA_index_smallest] )
      break;

    GOMEA_index_smallest++;
  }

  generationalStepAllGOMEAsRecursiveFold( GOMEA_index_smallest, GOMEA_index_biggest );
}

void generationalStepAllGOMEAsRecursiveFold( int GOMEA_index_smallest, int GOMEA_index_biggest )
{
  int i, GOMEA_index;

  for( i = 0; i < number_of_subgenerations_per_GOMEA_factor-1; i++ )
  {
    for( GOMEA_index = GOMEA_index_smallest; GOMEA_index <= GOMEA_index_biggest; GOMEA_index++ )
    {
      if( !terminated[GOMEA_index] )
      {
        terminated[GOMEA_index] = checkTerminationSpecificGOMEA( GOMEA_index );
      }

      if( (!terminated[GOMEA_index]) && (GOMEA_index >= minimum_GOMEA_index) )
      {
        makeOffspringSpecificGOMEA( GOMEA_index );

        selectFinalSurvivorsSpecificGOMEA( GOMEA_index );

        computeAverageFitnessSpecificGOMEA( GOMEA_index );

        if( earlystop )
          return;
      }
    }

    for( GOMEA_index = GOMEA_index_smallest; GOMEA_index < GOMEA_index_biggest; GOMEA_index++ )
      generationalStepAllGOMEAsRecursiveFold( GOMEA_index_smallest, GOMEA_index );
  }
}

void makeOffspringSpecificGOMEA( int gomea_index )
{
  // EDIT: RE-ENCODING OF SOLUTIONS
  if (reencode) {
    double *surrogate = (double *) Malloc( number_of_parameters*sizeof( double ) );
    for( int p = 0; p < population_sizes[gomea_index]; p++ )
    {
      for( int i = 0; i < randomkeysuntil; i++ )
        surrogate[i] = randomRealUniform01();
      int *order1 = mergeSortDoublesIncreasing( populations[gomea_index][p], randomkeysuntil );
      int *order2 = mergeSortDoublesIncreasing( surrogate, randomkeysuntil );
      for( int i = 0; i < randomkeysuntil; i++ )
        populations[gomea_index][p][order1[i]] = surrogate[order2[i]];
      free( order1 );
      free( order2 );
    }
    for( int i = 0; i < randomkeysuntil; i++ )
      surrogate[i] = randomRealUniform01();
    int *order1 = mergeSortDoublesIncreasing( elitist_solution, randomkeysuntil );
    int *order2 = mergeSortDoublesIncreasing( surrogate, randomkeysuntil );
    for( int i = 0; i < randomkeysuntil; i++ )
      elitist_solution[order1[i]] = surrogate[order2[i]];
    free( surrogate );
  }

  learnFOSSpecificGOMEA( gomea_index );

  if( print_FOSs_contents )
  {
    printf("### FOS contents for GOMEA #%02d in generation #%03d\n", gomea_index, number_of_generations);
    printFOSContentsSpecificGOMEA( gomea_index );
    printf( "###################################################\n" );
  }

  generateAndEvaluateNewSolutionsToFillOffspringSpecificGOMEA( gomea_index );
}

/**
 * Returns the name of an installed problem.
 */
char *installedFOSStructureName( int index )
{
  switch( index )
  {
    case  0: return( (char *) "Univariate" );
    case  1: return( (char *) "Linkage Tree" );
    case  2: return( (char *) "Multiscale Linkage Neighbors" );
    case  3: return( (char *) "Linkage Trees and Neighbors" );
    case  4: return( (char *) "Filtered Linkage Tree" );
    case  5: return( (char *) "Filtered Multiscale Linkage Neighbors" );
    case  6: return( (char *) "Filtered Linkage Trees and Neighbors" );
  }

  return( NULL );
}

/**
 * Returns the number of FOS structures installed.
 */
int numberOfInstalledFOSStructures()
{
  static int result = -1;
  
  if( result == -1 )
  {
    result = 0;
    while( installedFOSStructureName( result ) != NULL )
      result++;
  }
  
  return( result );
}

/**
 * Selects the FOS to be learned and calls the appropriate function to do
 * the learning.
 */
void learnFOSSpecificGOMEA( int gomea_index )
{
  int i;

  if( FOSs[gomea_index] != NULL )
  {
    for( i = 0; i < FOSs_length[gomea_index]; i++ )
      free( FOSs[gomea_index][i] );
    free( FOSs[gomea_index] );
    free( FOSs_number_of_indices[gomea_index] );
  }

  switch( FOSs_structure_index )
  {
    case  0: learnUnivariateFOSSpecificGOMEA( gomea_index ); break;
    case  1: learnLTFOSSpecificGOMEA( gomea_index, 1, 0, NULL ); break;
    case  2: learnMLNFOSSpecificGOMEA( gomea_index, 1, 0, NULL ); break;
    case  3: learnLTNFOSSpecificGOMEA( gomea_index ); break;
    case  4: learnFilteredLTFOSSpecificGOMEA( gomea_index, 1 ); break;
    case  5: learnFilteredMLNFOSSpecificGOMEA( gomea_index, 1 ); break;
    case  6: learnFilteredLTNFOSSpecificGOMEA( gomea_index ); break;
  }

  switch( FOSs_structure_index )
  {
    case  0: break;
    case  1: break;
    case  2: uniquifyFOSSpecificGOMEA( gomea_index ); break;
    case  3: uniquifyFOSSpecificGOMEA( gomea_index ); break;
    case  4: break;
    case  5: uniquifyFOSSpecificGOMEA( gomea_index ); break;
    case  6: uniquifyFOSSpecificGOMEA( gomea_index ); break;
    case  7: break;
  }
}

/**
 * Learns a univariate FOS (randomized ordering of the singletons).
 */
void learnUnivariateFOSSpecificGOMEA( int gomea_index )
{
  int i, FOSs_index, *order;

  order                               = randomPermutation( number_of_parameters );
  FOSs_length[gomea_index]            = number_of_parameters;
  FOSs[gomea_index]                   = (int **) Malloc( FOSs_length[gomea_index]*sizeof( int * ) );
  FOSs_number_of_indices[gomea_index] = (int *) Malloc( FOSs_length[gomea_index]*sizeof( int ) );
  FOSs_index                           = 0;
  for( i = 0; i < number_of_parameters; i++ )
  {
    FOSs[gomea_index][FOSs_index]                   = (int *) Malloc( 1*sizeof( int ) );
    FOSs[gomea_index][FOSs_index][0]                = order[FOSs_index];
    FOSs_number_of_indices[gomea_index][FOSs_index] = 1;
    FOSs_index++;
  }

  free( order );
}

/**
 * Learns a linkage tree FOS by means of hierarchical clustering.
 * This implementation follows the reciprocal nearest neighbor approach.
 */
int **learnLTFOSSpecificGOMEA( int gomea_index, short compute_dependency_matrices, short compute_parent_child_relations, int *number_of_parent_child_relations )
{
  char     done;
  int      i, j, r0, r1, rswap, *indices, *order,
           FOSs_index, **mpm, *mpm_number_of_indices, mpm_length,
         **mpm_new, *mpm_new_number_of_indices, mpm_new_length,
          *NN_chain, NN_chain_length, **parent_child_relations,
           PCR_index, *FOSs_index_of_mpm_element;
  double **S_matrix, mul0, mul1;

  parent_child_relations   = NULL; /* Only needed to prevent compiler warnings. */
  PCR_index                = 0;    /* Only needed to prevent compiler warnings. */
  FOSs_index_of_mpm_element = NULL; /* Only needed to prevent compiler warnings. */
  if( compute_parent_child_relations )
  {
    *number_of_parent_child_relations = number_of_parameters-1;
    parent_child_relations = (int **) Malloc( (*number_of_parent_child_relations)*sizeof( int * ) );
    for( i = 0; i < (*number_of_parent_child_relations); i++ )
      parent_child_relations[i] = (int *) Malloc( 3*sizeof( int ) );
    FOSs_index_of_mpm_element = (int *) Malloc( number_of_parameters*sizeof( int ) );
    for( i = 0; i < number_of_parameters; i++ )
      FOSs_index_of_mpm_element[i] = i;
  }

  /* Compute dependency matrix */
  if( compute_dependency_matrices )
    computeDependencyMatrixSpecificGOMEA( gomea_index );

  /* Initialize MPM to the univariate factorization */
  order                 = randomPermutation( number_of_parameters );
  mpm                   = (int **) Malloc( number_of_parameters*sizeof( int * ) );
  mpm_number_of_indices = (int *) Malloc( number_of_parameters*sizeof( int ) );
  mpm_length            = number_of_parameters;
  for( i = 0; i < number_of_parameters; i++ )
  {
    indices                  = (int *) Malloc( 1*sizeof( int ) );
    indices[0]               = order[i];
    mpm[i]                   = indices;
    mpm_number_of_indices[i] = 1;
  }
  free( order );

  /* Initialize LT to the initial MPM */
  FOSs_length[gomea_index]            = number_of_parameters+number_of_parameters-1;
  FOSs[gomea_index]                   = (int **) Malloc( FOSs_length[gomea_index]*sizeof( int * ) );
  FOSs_number_of_indices[gomea_index] = (int *) Malloc( FOSs_length[gomea_index]*sizeof( int ) );
  FOSs_index                                             = 0;
  for( i = 0; i < mpm_length; i++ )
  {
    FOSs[gomea_index][FOSs_index]                   = mpm[i];
    FOSs_number_of_indices[gomea_index][FOSs_index] = mpm_number_of_indices[i];
    FOSs_index++;
  }

  /* Initialize similarity matrix */
  S_matrix = (double **) Malloc( number_of_parameters*sizeof( double * ) );
  for( i = 0; i < number_of_parameters; i++ )
    S_matrix[i] = (double *) Malloc( number_of_parameters*sizeof( double ) );
  for( i = 0; i < mpm_length; i++ )
    for( j = 0; j < mpm_length; j++ )
      S_matrix[i][j] = dependency_matrices[gomea_index][mpm[i][0]][mpm[j][0]];
  for( i = 0; i < mpm_length; i++ )
    S_matrix[i][i] = 0;

  NN_chain        = (int *) Malloc( (number_of_parameters+2)*sizeof( int ) );
  NN_chain_length = 0;
  done            = 0;
  while( !done )
  {
    if( NN_chain_length == 0 )
    {
      NN_chain[NN_chain_length] = randomInt( mpm_length );
      NN_chain_length++;
    }

    while( NN_chain_length < 3 )
    {
      NN_chain[NN_chain_length] = determineNearestNeighbour( NN_chain[NN_chain_length-1], S_matrix, mpm_number_of_indices, mpm_length );
      NN_chain_length++;
    }

    while( NN_chain[NN_chain_length-3] != NN_chain[NN_chain_length-1] )
    {
      NN_chain[NN_chain_length] = determineNearestNeighbour( NN_chain[NN_chain_length-1], S_matrix, mpm_number_of_indices, mpm_length );
      if( ((S_matrix[NN_chain[NN_chain_length-1]][NN_chain[NN_chain_length]] == S_matrix[NN_chain[NN_chain_length-1]][NN_chain[NN_chain_length-2]])) && (NN_chain[NN_chain_length] != NN_chain[NN_chain_length-2]) )
        NN_chain[NN_chain_length] = NN_chain[NN_chain_length-2];
      NN_chain_length++;
      if( NN_chain_length > number_of_parameters )
        break;
    }
    r0 = NN_chain[NN_chain_length-2];
    r1 = NN_chain[NN_chain_length-1];
    if( r0 > r1 )
    {
      rswap = r0;
      r0    = r1;
      r1    = rswap;
    }
    NN_chain_length -= 3;

    if( r1 < mpm_length ) /* This test is required for exceptional cases in which the nearest-neighbor ordering has changed within the chain while merging within that chain */
    {
      indices = (int *) Malloc( (mpm_number_of_indices[r0]+mpm_number_of_indices[r1])*sizeof( int ) );
  
      i = 0;
      for( j = 0; j < mpm_number_of_indices[r0]; j++ )
      {
        indices[i] = mpm[r0][j];
        i++;
      }
      for( j = 0; j < mpm_number_of_indices[r1]; j++ )
      {
        indices[i] = mpm[r1][j];
        i++;
      }
    
      if( compute_parent_child_relations )
      {
        parent_child_relations[PCR_index][0] = FOSs_index;
        parent_child_relations[PCR_index][1] = FOSs_index_of_mpm_element[r0];
        parent_child_relations[PCR_index][2] = FOSs_index_of_mpm_element[r1];
        FOSs_index_of_mpm_element[r0]         = FOSs_index;
        FOSs_index_of_mpm_element[r1]         = FOSs_index_of_mpm_element[mpm_length-1];
        PCR_index++;
      }
      FOSs[gomea_index][FOSs_index]                   = indices;
      FOSs_number_of_indices[gomea_index][FOSs_index] = mpm_number_of_indices[r0]+mpm_number_of_indices[r1];
      FOSs_index++;
  
      mul0 = ((double) mpm_number_of_indices[r0])/((double) mpm_number_of_indices[r0]+mpm_number_of_indices[r1]);
      mul1 = ((double) mpm_number_of_indices[r1])/((double) mpm_number_of_indices[r0]+mpm_number_of_indices[r1]);
      for( i = 0; i < mpm_length; i++ )
      {
        if( (i != r0) && (i != r1) )
        {
          S_matrix[i][r0] = mul0*S_matrix[i][r0] + mul1*S_matrix[i][r1];
          S_matrix[r0][i] = S_matrix[i][r0];
        }
      }
  
      mpm_new                   = (int **) Malloc( (mpm_length-1)*sizeof( int * ) );
      mpm_new_number_of_indices = (int *) Malloc( (mpm_length-1)*sizeof( int ) );
      mpm_new_length            = mpm_length-1;
      for( i = 0; i < mpm_new_length; i++ )
      {
        mpm_new[i]                   = mpm[i];
        mpm_new_number_of_indices[i] = mpm_number_of_indices[i];
      }
  
      mpm_new[r0]                   = indices;
      mpm_new_number_of_indices[r0] = mpm_number_of_indices[r0]+mpm_number_of_indices[r1];
      if( r1 < mpm_length-1 )
      {
        mpm_new[r1]                   = mpm[mpm_length-1];
        mpm_new_number_of_indices[r1] = mpm_number_of_indices[mpm_length-1];
  
        for( i = 0; i < r1; i++ )
        {
          S_matrix[i][r1] = S_matrix[i][mpm_length-1];
          S_matrix[r1][i] = S_matrix[i][r1];
        }
  
        for( j = r1+1; j < mpm_new_length; j++ )
        {
          S_matrix[r1][j] = S_matrix[j][mpm_length-1];
          S_matrix[j][r1] = S_matrix[r1][j];
        }
      }
  
      for( i = 0; i < NN_chain_length; i++ )
      {
        if( NN_chain[i] == mpm_length-1 )
        {
          NN_chain[i] = r1;
          break;
        }
      }
  
      free( mpm );
      free( mpm_number_of_indices );
      mpm                   = mpm_new;
      mpm_number_of_indices = mpm_new_number_of_indices;
      mpm_length            = mpm_new_length;
  
      if( mpm_length == 1 )
        done = 1;
    }
  }

  free( NN_chain );

  free( mpm_new );
  free( mpm_number_of_indices );

  for( i = 0; i < number_of_parameters; i++ )
    free( S_matrix[i] );
  free( S_matrix );

  free( FOSs_index_of_mpm_element );

  return( parent_child_relations );
}

/**
 * Determines nearest neighbour according to similarity values.
 */
int determineNearestNeighbour( int index, double **S_matrix, int *mpm_number_of_indices, int mpm_length )
{
  int i, result;

  result = 0;
  if( result == index )
    result++;
  for( i = 1; i < mpm_length; i++ )
  {
    if( ((S_matrix[index][i] > S_matrix[index][result]) || ((S_matrix[index][i] == S_matrix[index][result]) && (mpm_number_of_indices[i] < mpm_number_of_indices[result]))) && (i != index) )
      result = i;
  }

  return( result );
}

/**
 * Learns a multiscale linkage neighbors FOS.
 */
int **learnMLNFOSSpecificGOMEA( int gomea_index, short compute_dependency_matrices, short compute_parent_child_relations, int *number_of_parent_child_relations )
{
  int    i, j, k, k2, k3, q, q2, **neighbors, ***buckets, **bucket_sizes, bucket_index, number_of_buckets,
         PCR_index, **parent_child_relations, **parent_child_relations_temp, number_of_parent_child_relations_temp, parent_child_relations_temp_size;
  double MI_max, MI_ratio;

  parent_child_relations           = NULL; /* Only needed to prevent compiler warnings. */
  parent_child_relations_temp      = NULL; /* Only needed to prevent compiler warnings. */
  PCR_index                        = 0;    /* Only needed to prevent compiler warnings. */
  parent_child_relations_temp_size = 0;    /* Only needed to prevent compiler warnings. */
  number_of_buckets                = 1+sqrt( population_sizes[gomea_index] );
  if( compute_parent_child_relations )
  {
    number_of_parent_child_relations_temp = 0;
    parent_child_relations_temp_size      = number_of_parameters*number_of_buckets;
    parent_child_relations_temp           = (int **) Malloc( (parent_child_relations_temp_size)*sizeof( int * ) );
    for( i = 0; i < parent_child_relations_temp_size; i++ )
      parent_child_relations_temp[i] = (int *) Malloc( 3*sizeof( int ) );
  }

  /* Compute dependency matrix */
  if( compute_dependency_matrices )
    computeDependencyMatrixSpecificGOMEA( gomea_index );

  /* Create a random ordering of neighbors for each variable */
  neighbors = (int **) Malloc( number_of_parameters*sizeof( int * ) );
  for( i = 0; i < number_of_parameters; i++ )
    neighbors[i] = randomPermutation( number_of_parameters );

  /* Determine for each variable i a particular ordering: a bucket sort without sorted buckets */
  buckets           = (int ***) Malloc( number_of_parameters*sizeof( int ** ) );
  bucket_sizes      = (int **) Malloc( number_of_parameters*sizeof( int * ) );
  for( i = 0; i < number_of_parameters; i++ )
  {
    buckets[i]      = NULL;
    bucket_sizes[i] = NULL;

    buckets[i]      = (int **) Malloc( number_of_buckets*sizeof( int * ) );
    bucket_sizes[i] = (int *) Malloc( number_of_buckets*sizeof( int ) );

    MI_max = 0;
    for( j = 0; j < number_of_parameters; j++ )
    {
      if( neighbors[i][j] != i )
      {
        if( dependency_matrices[gomea_index][i][neighbors[i][j]] > MI_max )
          MI_max = dependency_matrices[gomea_index][i][neighbors[i][j]];
      }
    }
    MI_max = MI_max > 1 ? 1 : MI_max;

    for( bucket_index = 0; bucket_index < number_of_buckets; bucket_index++ )
      bucket_sizes[i][bucket_index] = 0;
    for( j = 0; j < number_of_parameters; j++ )
    {
      if( neighbors[i][j] == i )
        bucket_index = number_of_buckets-1;
      else
      {
        MI_ratio = 1;
        if( MI_max > 0 )
          MI_ratio = dependency_matrices[gomea_index][i][neighbors[i][j]]/MI_max;
        bucket_index = MI_ratio*((double) (number_of_buckets-1));
        bucket_index = bucket_index >= (number_of_buckets-1) ? (number_of_buckets-1)-1 : bucket_index;
        bucket_index = bucket_index <= 0 ? 0 : bucket_index;
      }
      bucket_sizes[i][bucket_index]++;
    }

    for( bucket_index = 0; bucket_index < number_of_buckets; bucket_index++ )
    {
      buckets[i][bucket_index] = NULL;
      if( bucket_sizes[i][bucket_index] > 0 )
        buckets[i][bucket_index] = (int *) Malloc( bucket_sizes[i][bucket_index]*sizeof( int ) );
    }

    for( bucket_index = 0; bucket_index < number_of_buckets; bucket_index++ )
      bucket_sizes[i][bucket_index] = 0;
    for( j = 0; j < number_of_parameters; j++ )
    {
      if( neighbors[i][j] == i )
        bucket_index = number_of_buckets-1;
      else
      {
        MI_ratio = 1;
        if( MI_max > 0 )
          MI_ratio = dependency_matrices[gomea_index][i][neighbors[i][j]]/MI_max;
        bucket_index = MI_ratio*((double) (number_of_buckets-1));
        bucket_index = bucket_index >= (number_of_buckets-1) ? (number_of_buckets-1)-1 : bucket_index;
        bucket_index = bucket_index <= 0 ? 0 : bucket_index;
      }
      buckets[i][bucket_index][bucket_sizes[i][bucket_index]] = neighbors[i][j];
      bucket_sizes[i][bucket_index]++;
    }
  }

  FOSs_length[gomea_index] = 0;
  for( i = 0; i < number_of_parameters; i++ )
  {
    if( number_of_parameters == 1 )
      FOSs_length[gomea_index]++;
    else
    {
      if( bucket_sizes[i][number_of_buckets-1] > 1 )
        FOSs_length[gomea_index]++;
      q  = bucket_sizes[i][number_of_buckets-1];
      q2 = q;
      for( j = number_of_buckets-2; j >= 0; j-- )
      {
        q2 += bucket_sizes[i][j];
        if( ((q2 >= (2*q)) || (q2 == number_of_parameters)) && (q2 <= (number_of_parameters/2)) )
        {
          q = q2;
          if( compute_parent_child_relations )
          {
            parent_child_relations_temp[PCR_index][0] = FOSs_length[gomea_index];
            parent_child_relations_temp[PCR_index][1] = FOSs_length[gomea_index]-1;
            parent_child_relations_temp[PCR_index][2] = -1;
            PCR_index++;
          }
          FOSs_length[gomea_index]++;
        }
        if( q2 == number_of_parameters )
          break;
      }
    }
  }

  /* Create the multiscale set of neighbors for each variable */
  FOSs[gomea_index]                   = (int **) Malloc( FOSs_length[gomea_index]*sizeof( int * ) );
  FOSs_number_of_indices[gomea_index] = (int *) Malloc( FOSs_length[gomea_index]*sizeof( int ) );
  
  FOSs_length[gomea_index] = 0;
  for( i = 0; i < number_of_parameters; i++ )
  {
    if( number_of_parameters == 1 )
    {
      FOSs_number_of_indices[gomea_index][FOSs_length[gomea_index]] = 1;
      FOSs[gomea_index][FOSs_length[gomea_index]]                   = (int *) Malloc( FOSs_number_of_indices[gomea_index][FOSs_length[gomea_index]]*sizeof( int ) );
      FOSs[gomea_index][FOSs_length[gomea_index]][0]                = neighbors[i][0];
      FOSs_length[gomea_index]++;
    }
    else
    {
      if( bucket_sizes[i][number_of_buckets-1] > 1 )
      {
        FOSs_number_of_indices[gomea_index][FOSs_length[gomea_index]] = bucket_sizes[i][number_of_buckets-1];
        FOSs[gomea_index][FOSs_length[gomea_index]]                   = (int *) Malloc( FOSs_number_of_indices[gomea_index][FOSs_length[gomea_index]]*sizeof( int ) );
        for( k = 0; k < bucket_sizes[i][number_of_buckets-1]; k++ )
          FOSs[gomea_index][FOSs_length[gomea_index]][k] = buckets[i][number_of_buckets-1][k];
        FOSs_length[gomea_index]++;
      }
      q  = bucket_sizes[i][number_of_buckets-1];
      q2 = q;
      for( j = number_of_buckets-2; j >= 0; j-- )
      {
        q2 += bucket_sizes[i][j];
        if( ((q2 >= (2*q)) || (q2 == number_of_parameters)) && (q2 <= (number_of_parameters/2)) )
        {
          q = q2;
          FOSs_number_of_indices[gomea_index][FOSs_length[gomea_index]] = q;
          FOSs[gomea_index][FOSs_length[gomea_index]]                   = (int *) Malloc( FOSs_number_of_indices[gomea_index][FOSs_length[gomea_index]]*sizeof( int ) );
          k                                 = 0;
          for( k2 = number_of_buckets-1; k2 >= j; k2-- )
          {
            for( k3 = 0; k3 < bucket_sizes[i][k2]; k3++ )
            {
              FOSs[gomea_index][FOSs_length[gomea_index]][k] = buckets[i][k2][k3];
              k++;
            }
          }
          FOSs_length[gomea_index]++;
        }
        if( q2 == number_of_parameters )
          break;
      }
    }
  }

  for( i = 0; i < number_of_parameters; i++ )
  {
    if( buckets[i] != NULL )
    {
      for( j = 0; j < number_of_buckets; j++ )
      {
        if( buckets[i][j] != NULL )
        {
          free( buckets[i][j] );
        }
      }
      free( buckets[i] );
      free( bucket_sizes[i] );
    }
  }
  free( bucket_sizes );
  free( buckets );

  for( i = 0; i < number_of_parameters; i++ )
  {
    if( neighbors[i] != NULL )
    {
      free( neighbors[i] );
    }
  }
  free( neighbors );

  if( compute_parent_child_relations )
  {
    number_of_parent_child_relations_temp = PCR_index;
    *number_of_parent_child_relations     = number_of_parent_child_relations_temp;

    parent_child_relations = (int **) Malloc( (*number_of_parent_child_relations)*sizeof( int * ) );
    for( i = 0; i < (*number_of_parent_child_relations); i++ )
    {
      parent_child_relations[i] = (int *) Malloc( 3*sizeof( int ) );
      for( j = 0; j < 3; j++ )
        parent_child_relations[i][j] = parent_child_relations_temp[i][j];
    }
    for( i = 0; i < parent_child_relations_temp_size; i++ )
      free( parent_child_relations_temp[i] );
    free( parent_child_relations_temp );
  }
  return( parent_child_relations );
}

/**
 * Learns a multiscale linkage neighbors FOS.
 */
void learnLTNFOSSpecificGOMEA( int gomea_index )
{
  learnLTNFOSWithOrWithoutFilteringSpecificGOMEA( gomea_index, 0 );
}

/**
 * Learns a multiscale linkage neighbors FOS with or without filtering. The actual
 * work is done here.
 */
void learnLTNFOSWithOrWithoutFilteringSpecificGOMEA( int gomea_index, short use_filtering )
{
  int i, j, **LT_FOS, **MLN_FOS, *LT_FOS_number_of_indices, *MLN_FOS_number_of_indices, LT_FOS_length, MLN_FOS_length;

  /* Learn the LT FOS and create a backup copy */
  if( use_filtering )
    learnFilteredLTFOSSpecificGOMEA( gomea_index, 1 );
  else
    learnLTFOSSpecificGOMEA( gomea_index, 1, 0, NULL );
  
  LT_FOS_length            = FOSs_length[gomea_index];
  LT_FOS_number_of_indices = (int *) Malloc( FOSs_length[gomea_index]*sizeof( int ) );
  LT_FOS                   = (int **) Malloc( FOSs_length[gomea_index]*sizeof( int * ) );
  for( i = 0; i < FOSs_length[gomea_index]; i++ )
  {
    LT_FOS_number_of_indices[i] = FOSs_number_of_indices[gomea_index][i];
    LT_FOS[i] = (int *) Malloc( FOSs_number_of_indices[gomea_index][i]*sizeof( int ) );
    for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
      LT_FOS[i][j] = FOSs[gomea_index][i][j];
  }

  for( i = 0; i < FOSs_length[gomea_index]; i++ )
    free( FOSs[gomea_index][i] );
  free( FOSs[gomea_index] );
  free( FOSs_number_of_indices[gomea_index] );

  /* Learn the MLN FOS and create a backup copy */
  if( use_filtering )
    learnFilteredMLNFOSSpecificGOMEA( gomea_index, 0 );
  else
    learnMLNFOSSpecificGOMEA( gomea_index, 0, 0, NULL );

  MLN_FOS_length            = FOSs_length[gomea_index];
  MLN_FOS_number_of_indices = (int *) Malloc( FOSs_length[gomea_index]*sizeof( int ) );
  MLN_FOS                   = (int **) Malloc( FOSs_length[gomea_index]*sizeof( int * ) );
  for( i = 0; i < FOSs_length[gomea_index]; i++ )
  {
    MLN_FOS_number_of_indices[i] = FOSs_number_of_indices[gomea_index][i];
    MLN_FOS[i] = (int *) Malloc( FOSs_number_of_indices[gomea_index][i]*sizeof( int ) );
    for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
      MLN_FOS[i][j] = FOSs[gomea_index][i][j];
  }

  for( i = 0; i < FOSs_length[gomea_index]; i++ )
    free( FOSs[gomea_index][i] );
  free( FOSs[gomea_index] );
  free( FOSs_number_of_indices[gomea_index] );

  /* Construct the LTN FOS: join the LT FOS and the MLN FOS */
  FOSs_length[gomea_index]            = LT_FOS_length + MLN_FOS_length; 
  FOSs_number_of_indices[gomea_index] = (int *) Malloc( (LT_FOS_length + MLN_FOS_length)*sizeof( int ) );
  FOSs[gomea_index]                   = (int **) Malloc( (LT_FOS_length + MLN_FOS_length)*sizeof( int * ) );
  for( i = 0; i < LT_FOS_length; i++ )
  {
    FOSs_number_of_indices[gomea_index][i] = LT_FOS_number_of_indices[i];
    FOSs[gomea_index][i] = (int *) Malloc( LT_FOS_number_of_indices[i]*sizeof( int ) );
    for( j = 0; j < LT_FOS_number_of_indices[i]; j++ )
      FOSs[gomea_index][i][j] = LT_FOS[i][j];
  }
  for( i = 0; i < MLN_FOS_length; i++ )
  {
    FOSs_number_of_indices[gomea_index][LT_FOS_length+i] = MLN_FOS_number_of_indices[i];
    FOSs[gomea_index][LT_FOS_length+i] = (int *) Malloc( MLN_FOS_number_of_indices[i]*sizeof( int ) );
    for( j = 0; j < MLN_FOS_number_of_indices[i]; j++ )
      FOSs[gomea_index][LT_FOS_length+i][j] = MLN_FOS[i][j];
  }

  /* Free up backup memory */
  for( i = 0; i < LT_FOS_length; i++ )
    free( LT_FOS[i] );
  free( LT_FOS_number_of_indices );
  free( LT_FOS );
  for( i = 0; i < MLN_FOS_length; i++ )
    free( MLN_FOS[i] );
  free( MLN_FOS_number_of_indices );
  free( MLN_FOS );
}

void learnFilteredLTFOSSpecificGOMEA( int gomea_index, short compute_dependency_matrices )
{
  int i, **parent_child_relations, number_of_parent_child_relations;

  parent_child_relations = learnLTFOSSpecificGOMEA( gomea_index, compute_dependency_matrices, 1, &number_of_parent_child_relations );

  filterParentChildRelationsAndCreateNewFOSSpecificGOMEA( gomea_index, parent_child_relations, number_of_parent_child_relations );

  for( i = 0; i < number_of_parent_child_relations; i++ )
    free( parent_child_relations[i] );
  free( parent_child_relations );
}

void learnFilteredMLNFOSSpecificGOMEA( int gomea_index, short compute_dependency_matrices )
{
  int i, **parent_child_relations, number_of_parent_child_relations;

  parent_child_relations = learnMLNFOSSpecificGOMEA( gomea_index, compute_dependency_matrices, 1, &number_of_parent_child_relations );

  filterParentChildRelationsAndCreateNewFOSSpecificGOMEA( gomea_index, parent_child_relations, number_of_parent_child_relations );

  for( i = 0; i < number_of_parent_child_relations; i++ )
    free( parent_child_relations[i] );
  free( parent_child_relations );
}

/**
 * Learns a multiscale linkage neighbors FOS.
 */
void learnFilteredLTNFOSSpecificGOMEA( int gomea_index )
{
  learnLTNFOSWithOrWithoutFilteringSpecificGOMEA( gomea_index, 1 );
}

void filterParentChildRelationsAndCreateNewFOSSpecificGOMEA( int gomea_index, int **parent_child_relations, int number_of_parent_child_relations )
{
  char   *FOS_element_accepted, *linkage_strength_computed;
  int     i, j, parent_index, child0_index, child1_index, **FOS, *FOS_number_of_indices, FOS_length;
  double *linkage_strength;

  FOS_element_accepted = (char *) Malloc( FOSs_length[gomea_index]*sizeof( char ) );
  for( i = 0; i < FOSs_length[gomea_index]; i++ )
    FOS_element_accepted[i] = 1;

  linkage_strength_computed = (char *) Malloc( FOSs_length[gomea_index]*sizeof( char ) );
  for( i = 0; i < FOSs_length[gomea_index]; i++ )
    linkage_strength_computed[i] = 0;

  linkage_strength = (double *) Malloc( FOSs_length[gomea_index]*sizeof( double ) );
  for( i = 0; i < FOSs_length[gomea_index]; i++ )
    linkage_strength[i] = 0;

  for( i = 0; i < number_of_parent_child_relations; i++ )
  {
    parent_index = parent_child_relations[i][0];
    child0_index = parent_child_relations[i][1];
    child1_index = parent_child_relations[i][2];

    if( !linkage_strength_computed[parent_index] )
    {
      linkage_strength[parent_index]          = computeLinkageStrengthSpecificGOMEA( gomea_index, FOSs[gomea_index][parent_index], FOSs_number_of_indices[gomea_index][parent_index] );
      linkage_strength_computed[parent_index] = 1;
    }
    if( child0_index != -1 )
    {
      if( !linkage_strength_computed[child0_index] )
      {
        linkage_strength[child0_index]          = computeLinkageStrengthSpecificGOMEA( gomea_index, FOSs[gomea_index][child0_index], FOSs_number_of_indices[gomea_index][child0_index] );
        linkage_strength_computed[child0_index] = 1;
      }
    }
    if( child1_index != -1 )
    {
      if( !linkage_strength_computed[child1_index] )
      {
        linkage_strength[child1_index]          = computeLinkageStrengthSpecificGOMEA( gomea_index, FOSs[gomea_index][child1_index], FOSs_number_of_indices[gomea_index][child1_index] );
        linkage_strength_computed[child1_index] = 1;
      }
    }

    /* Remove each child if it has the same linkage strength as its parent */
    if( child0_index != -1 )
    {
      if( ((linkage_strength[parent_index] >= (1.0-1e-100)*linkage_strength[child0_index])) && ((linkage_strength[parent_index] <= (1.0+1e-100)*linkage_strength[child0_index])) )
        FOS_element_accepted[child0_index] = 0;
    }
    if( child1_index != -1 )
    {
      if( ((linkage_strength[parent_index] >= (1.0-1e-100)*linkage_strength[child1_index])) && ((linkage_strength[parent_index] <= (1.0+1e-100)*linkage_strength[child1_index])) )
        FOS_element_accepted[child1_index] = 0;
    }
  }

  /* Create a backup copy of the FOS */
  FOS_length            = FOSs_length[gomea_index];
  FOS_number_of_indices = (int *) Malloc( FOSs_length[gomea_index]*sizeof( int ) );
  FOS                   = (int **) Malloc( FOSs_length[gomea_index]*sizeof( int * ) );
  for( i = 0; i < FOSs_length[gomea_index]; i++ )
  {
    FOS_number_of_indices[i] = FOSs_number_of_indices[gomea_index][i];
    FOS[i] = (int *) Malloc( FOSs_number_of_indices[gomea_index][i]*sizeof( int ) );
    for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
      FOS[i][j] = FOSs[gomea_index][i][j];
  }

  for( i = 0; i < FOSs_length[gomea_index]; i++ )
    free( FOSs[gomea_index][i] );
  free( FOSs[gomea_index] );
  free( FOSs_number_of_indices[gomea_index] );

  /* Apply filters */
  FOSs_length[gomea_index] = 0;
  for( i = 0; i < FOS_length; i++ )
    if( FOS_element_accepted[i] )
      FOSs_length[gomea_index]++;
 
  FOSs_number_of_indices[gomea_index] = (int *) Malloc( FOSs_length[gomea_index]*sizeof( int ) );
  FOSs[gomea_index]                   = (int **) Malloc( FOSs_length[gomea_index]*sizeof( int * ) );
  FOSs_length[gomea_index]            = 0;
  for( i = 0; i < FOS_length; i++ )
  {
    if( FOS_element_accepted[i] )
    {
      FOSs_number_of_indices[gomea_index][FOSs_length[gomea_index]] = FOS_number_of_indices[i];
      FOSs[gomea_index][FOSs_length[gomea_index]]                   = (int *) Malloc( FOS_number_of_indices[i]*sizeof( int ) );
      for( j = 0; j < FOS_number_of_indices[i]; j++ )
        FOSs[gomea_index][FOSs_length[gomea_index]][j] = FOS[i][j];
      FOSs_length[gomea_index]++;
    }
  }

  /* Free up backup memory */
  for( i = 0; i < FOS_length; i++ )
    free( FOS[i] );
  free( FOS_number_of_indices );
  free( FOS );

  free( linkage_strength );
  free( linkage_strength_computed );
  free( FOS_element_accepted );
}

double computeLinkageStrengthSpecificGOMEA( int gomea_index, int *variables, int number_of_variables )
{
  int    i, j, n;
  double result;

  result = dependency_matrices[gomea_index][variables[0]][variables[0]];
  if( number_of_variables > 1 )
  {
    result = 0;
    n      = 0;
    for( i = 0; i < number_of_variables; i++ )
    {
      for( j = i+1; j < number_of_variables; j++ )
      {
        result += dependency_matrices[gomea_index][variables[i]][variables[j]];
        n++;
      }
    }
    result /= (double) n;
  }

  return( result );
}

void computeDependencyMatrixSpecificGOMEA( int gomea_index )
{
  int    i, j, k;
  double p, entropy, average_distance;

  if(dependency_matrix_method == 1 || dependency_matrix_method == 3) {
    // Keydistance
    // NOTE: PROBABLY ONLY WORKS WELL IF YOU CAN 'TRANSLATE' blocks.
    for( i = 0; i < number_of_parameters; i++ )
    {
      for( j = i+1; j < number_of_parameters; j++ )
      {
        /* Multiply by inverted average distance between variables */
        average_distance = 0.0;
        for( k = 0; k < population_sizes[gomea_index]; k++ )
        {
          average_distance += abs(populations[gomea_index][k][i] - populations[gomea_index][k][j]); // *(populations[gomea_index][k][i] - populations[gomea_index][k][j]);
        }
        average_distance /= (double) population_sizes[gomea_index];

        if (dependency_matrix_method == 1) {
          dependency_matrices[gomea_index][i][j] = average_distance;
        } else {
          dependency_matrices[gomea_index][i][j] = 1 - average_distance;
        }

        /* Create symmetric matrix */
        dependency_matrices[gomea_index][j][i] = dependency_matrices[gomea_index][i][j];
      }
    }
  } else if (dependency_matrix_method == 2) {
    for( i = 0; i < number_of_parameters; i++ )
    {
      for( j = i+1; j < number_of_parameters; j++ )
      {
      // Random dependency matrix.
      dependency_matrices[gomea_index][i][j] = randomRealUniform01();

      /* Create symmetric matrix */
      dependency_matrices[gomea_index][j][i] = dependency_matrices[gomea_index][i][j];
      }
    }
  } else {
    // EDIT: Original approach.
    for( i = 0; i < number_of_parameters; i++ )
    {
      for( j = i+1; j < number_of_parameters; j++ )
      {
        /* Compute entropy of probability that variable i has a smaller value than variable j, use inverted entropy */
        p = 0.0;
        for( k = 0; k < population_sizes[gomea_index]; k++ )
        {
          if( populations[gomea_index][k][i] < populations[gomea_index][k][j] )
            p += 1.0;
        }
        p       /= (double) population_sizes[gomea_index];
        entropy  = (p==0)||(p==1)?0:-(p*log2(p) + (1.0-p)*log2(1.0-p));

        dependency_matrices[gomea_index][i][j] = 1.0-entropy;

        /* Multiply by inverted average distance between variables */
        average_distance = 0.0;
        for( k = 0; k < population_sizes[gomea_index]; k++ )
        {
          average_distance += (populations[gomea_index][k][i] - populations[gomea_index][k][j])*(populations[gomea_index][k][i] - populations[gomea_index][k][j]);
        }
        average_distance /= (double) population_sizes[gomea_index];

        dependency_matrices[gomea_index][i][j] *= 1.0-average_distance;

        /* Create symmetric matrix */
        dependency_matrices[gomea_index][j][i] = dependency_matrices[gomea_index][i][j];
      }
    }
  }
}

/**
 * Estimates the cumulative probability distribution of a
 * single binary marginal.
 */
double *estimateParametersForSingleBinaryMarginalSpecificGOMEA( int gomea_index, int *indices, int number_of_indices, int *factor_size )
{
  int     i, j, index, power_of_two;
  double *result;

  *factor_size = (int) pow( 2, number_of_indices );
  result       = (double *) Malloc( (*factor_size)*sizeof( double ) );

  for( i = 0; i < (*factor_size); i++ )
    result[i] = 0.0;

  for( i = 0; i < population_sizes[gomea_index]; i++ )
  {
    index        = 0;
    power_of_two = 1;
    for( j = number_of_indices-1; j >= 0; j-- )
    {
      index += populations[gomea_index][i][indices[j]] ? power_of_two : 0;
      power_of_two *= 2;
    }

    result[index] += 1.0;
  }

  for( i = 0; i < (*factor_size); i++ )
    result[i] /= (double) population_sizes[gomea_index];

  for( i = 1; i < (*factor_size); i++ )
    result[i] += result[i-1];

  result[(*factor_size)-1] = 1.0;

  return( result );
}

/**
 * Ensures that every FOS element is unique, i.e. there are no
 * duplicate linkage subsets.
 */
void uniquifyFOSSpecificGOMEA( int gomea_index )
{
  short *FOSs_subset_is_duplicate;
  int    i, j, k, q, *sorted, **FOSs_new, *FOSs_number_of_indices_new, FOSs_length_new;

  /* First analyze which subsets are duplicates */
  FOSs_subset_is_duplicate = (short *) Malloc( FOSs_length[gomea_index]*sizeof( short ) );
  for( i = 0; i < FOSs_length[gomea_index]; i++ )
    FOSs_subset_is_duplicate[i] = 0;

  sorted = mergeSortIntegersDecreasing( FOSs_number_of_indices[gomea_index], FOSs_length[gomea_index] );

  i = 0;
  while( i < FOSs_length[gomea_index] )
  {
    /* Determine stretch of FOS elements that have the same length */
    j = i+1;
    while( (j < FOSs_length[gomea_index]) && (FOSs_number_of_indices[gomea_index][sorted[i]] == FOSs_number_of_indices[gomea_index][sorted[j]]) )
      j++;

    /* Check inside stretch for duplicates */
    for( k = i; k < j-1; k++ )
    {
      if( FOSs_subset_is_duplicate[sorted[k]] )
        continue;

      for( q = k+1; q < j; q++ )
      {
        if( FOSs_subset_is_duplicate[sorted[q]] )
          continue;

        if( linkageSubsetsOfSameLengthAreDuplicates( FOSs[gomea_index][sorted[k]], FOSs[gomea_index][sorted[q]], FOSs_number_of_indices[gomea_index][sorted[k]] ) )
          FOSs_subset_is_duplicate[sorted[q]] = 1;
      }
    }
    i = j;
  }

  /* Then re-create the FOS without the duplicate sets */
  FOSs_length_new = 0;
  for( i = 0; i < FOSs_length[gomea_index]; i++ )
  {
    if( !FOSs_subset_is_duplicate[i] )
      FOSs_length_new++;
  }

  FOSs_new                   = (int **) Malloc( FOSs_length_new*sizeof( int * ) );
  FOSs_number_of_indices_new = (int *) Malloc( FOSs_length_new*sizeof( int ) );

  j = 0;
  for( i = 0; i < FOSs_length[gomea_index]; i++ )
  {
    if( !FOSs_subset_is_duplicate[i] )
    {
      FOSs_new[j] = (int *) Malloc( FOSs_number_of_indices[gomea_index][i]*sizeof( int ) );
      for( k = 0; k < FOSs_number_of_indices[gomea_index][i]; k++ )
        FOSs_new[j][k] = FOSs[gomea_index][i][k];

      FOSs_number_of_indices_new[j] = FOSs_number_of_indices[gomea_index][i];

      j++;
    }
  }

  for( i = 0; i < FOSs_length[gomea_index]; i++ )
    free( FOSs[gomea_index][i] );
  free( FOSs[gomea_index] );
  free( FOSs_number_of_indices[gomea_index] );
  FOSs[gomea_index]                   = FOSs_new;
  FOSs_number_of_indices[gomea_index] = FOSs_number_of_indices_new;
  FOSs_length[gomea_index]            = FOSs_length_new;

  free( sorted );
  free( FOSs_subset_is_duplicate );
}

short linkageSubsetsOfSameLengthAreDuplicates( int *linkageSubset0, int *linkageSubset1, int length )
{
  short result, *linkageSubset0AsBitString, *linkageSubset1AsBitString;
  int   i, *sorted0, *sorted1;

  result = 0;
  if( length == 1 )
  {
    if( linkageSubset0[0] == linkageSubset1[0] )
      result = 1;
  }
  else if( length == number_of_parameters )
  {
    result = 1;
  }
  else if( length <= (number_of_parameters/2) )
  {
    sorted0 = mergeSortIntegersDecreasing( linkageSubset0, length );
    sorted1 = mergeSortIntegersDecreasing( linkageSubset1, length );
    result = 1;
    for( i = 0; i < length; i++ )
    {
      if( linkageSubset0[sorted0[i]] != linkageSubset1[sorted1[i]] )
      {
        result = 0;
        break;
      }
    }
    free( sorted0 );
    free( sorted1 );
  }
  else
  {
    linkageSubset0AsBitString = (short *) Malloc( number_of_parameters*sizeof( short ) );
    linkageSubset1AsBitString = (short *) Malloc( number_of_parameters*sizeof( short ) );
    for( i = 0; i < number_of_parameters; i++ )
    {
      linkageSubset0AsBitString[i] = 0;
      linkageSubset1AsBitString[i] = 0;
    }
    for( i = 0; i < length; i++ )
    {
      linkageSubset0AsBitString[linkageSubset0[i]] = 1;
      linkageSubset1AsBitString[linkageSubset1[i]] = 1;
    }
    result = 1;
    for( i = 0; i < number_of_parameters; i++ )
    {
      if( linkageSubset0AsBitString[i] != linkageSubset1AsBitString[i] )
      {
        result = 0;
        break;
      }
    }
    free( linkageSubset0AsBitString );
    free( linkageSubset1AsBitString );
  }

  return( result );
}

void printFOSContentsSpecificGOMEA( int gomea_index )
{
  int i, j;

  for( i = 0; i < FOSs_length[gomea_index]; i++ )
  {
    printf( "# [" );
    for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
    {
      printf( "%d",FOSs[gomea_index][i][j] );
      if( j < FOSs_number_of_indices[gomea_index][i]-1 )
        printf( " " );
    }
    printf( "]\n" );
  }
  fflush( stdout );
}

/**
 * Computes the two-log of x.
 */
double math_log_two = log(2.0);
double log2( double x )
{
  return( log(x) / math_log_two );
}

/**
 * Generates new solutions.
 */
void generateAndEvaluateNewSolutionsToFillOffspringSpecificGOMEA( int gomea_index )
{
  int    i, j;
  double obj, con, *solution;

  for( i = 0; i < population_sizes[gomea_index]; i++ )
  {
    solution = generateAndEvaluateNewSolutionSpecificGOMEA( gomea_index, i%(population_sizes[gomea_index]), &obj, &con );

    for( j = 0; j < number_of_parameters; j++ )
      offsprings[gomea_index][i][j] = solution[j];

    objective_values_offsprings[gomea_index][i] = obj;
    constraint_values_offsprings[gomea_index][i] = con;

    free( solution );

    if( earlystop )
      return;
  }
}

/**
 * Performs Genepool Optimal Mixing (for one solution in the population).
 */
double *generateAndEvaluateNewSolutionSpecificGOMEA( int gomea_index, int parent_index, double *obj, double *con )
{
  char    donor_parameters_are_the_same;
  double *result, *backup, *backup_for_ilse;
  short   solution_has_changed;
  int     i, j, j_ilse_best, index, parameter_index_for_ilse[1];
  double  obj_backup, con_backup, obj_backup_for_ilse, con_backup_for_ilse, obj_ilse_best, con_ilse_best;

  solution_has_changed = 0;

  result = (double *) Malloc( number_of_parameters*sizeof( double ) );
  for( i = 0; i < number_of_parameters; i++ )
    result[i] = populations[gomea_index][parent_index][i];

  *obj = objective_values[gomea_index][parent_index];
  *con = constraint_values[gomea_index][parent_index];

  backup = (double *) Malloc( number_of_parameters*sizeof( double ) );
  for( i = 0; i < number_of_parameters; i++ )
    backup[i] = result[i];

  backup_for_ilse = NULL; /* Only needed to prevent compiler warnings. */
  if( use_ilse )
  {
    backup_for_ilse = (double *) Malloc( number_of_parameters*sizeof( double ) );
    for( i = 0; i < number_of_parameters; i++ )
      backup_for_ilse[i] = result[i];
  }

  obj_backup = *obj;
  con_backup = *con;

  /* Phase 1: optimal mixing with random donors */
  shuffleFOSSpecificGOMEA( gomea_index );
  if( use_ilse )
    shuffleFOSSubsetsSpecificGOMEA( gomea_index );

  for( i = 0; i < FOSs_length[gomea_index]; i++ )
  {
    if( FOSs_number_of_indices[gomea_index][i] == number_of_parameters )
      continue;

    index = randomInt( population_sizes[gomea_index] );

    /* Convert index to binary representation and set factor variables. */
    for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
      result[FOSs[gomea_index][i][j]] = populations[gomea_index][index][FOSs[gomea_index][i][j]];

    // EDIT: RANDOM RESCALING
    if( rescale && randomRealUniform01() < 0.1)
    {
      int    number_of_intervals;
      double min, max, range, random_left_bracket;
      short  anykeys = 0;

      number_of_intervals = number_of_parameters;

      min = result[FOSs[gomea_index][i][0]]; 
      max = result[FOSs[gomea_index][i][0]]; 
      
      for( j = 1; j < FOSs_number_of_indices[gomea_index][i]; j++ )
      {
        // EDIT: Do not take into account non-random keys for rescaling.
        if( FOSs[gomea_index][i][j] >= randomkeysuntil )
          continue;
        if( result[FOSs[gomea_index][i][j]] < min )
          min = result[FOSs[gomea_index][i][j]];
        if( result[FOSs[gomea_index][i][j]] > max )
          max = result[FOSs[gomea_index][i][j]];
        anykeys = 1;
      }
      if( anykeys ) {
        range = max-min;

        random_left_bracket = ((double) randomInt(number_of_intervals))/((double)number_of_intervals);
        if( range > 0 )
        {
          for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
            if( FOSs[gomea_index][i][j] < randomkeysuntil )
              result[FOSs[gomea_index][i][j]] = ((result[FOSs[gomea_index][i][j]]-min)/range)*(1.0/((double)number_of_intervals)) + random_left_bracket;
        }
        else
        {
          for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
            if( FOSs[gomea_index][i][j] < randomkeysuntil )
              result[FOSs[gomea_index][i][j]] = random_left_bracket;
        }
      }
    }
    if(translate && randomRealUniform01() < 0.1) 
    {
      double delta;
      int v;
      v = randomInt(FOSs_number_of_indices[gomea_index][i]);

      delta = (2.0 * randomRealUniform01() - 0.5) * ((backup[FOSs[gomea_index][i][v]]) - (result[FOSs[gomea_index][i][v]]));

      for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ ) {
        result[FOSs[gomea_index][i][j]] += delta;
      }
    }
    
    /* Test if the change is for the better */
    donor_parameters_are_the_same = 1;
    for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
    {
      if( backup[FOSs[gomea_index][i][j]] != result[FOSs[gomea_index][i][j]] )
      {
        donor_parameters_are_the_same = 0;
        break;
      }
    }
    if( !donor_parameters_are_the_same )
    {
      if( !use_ilse )
        installedProblemEvaluation( problem_index, result, obj, con, FOSs_number_of_indices[gomea_index][i], FOSs[gomea_index][i], backup, obj_backup, con_backup );
      else
      {
        for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
          result[FOSs[gomea_index][i][j]] = backup[FOSs[gomea_index][i][j]];

        j_ilse_best         = 0;
        obj_ilse_best       = *obj;
        con_ilse_best       = *con;
        obj_backup_for_ilse = *obj;
        con_backup_for_ilse = *con;
        for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
        {
          if( result[FOSs[gomea_index][i][j]] != populations[gomea_index][index][FOSs[gomea_index][i][j]] )
          {
            result[FOSs[gomea_index][i][j]] = populations[gomea_index][index][FOSs[gomea_index][i][j]];
            parameter_index_for_ilse[0] = FOSs[gomea_index][i][j];
            installedProblemEvaluation( problem_index, result, obj, con, 1, parameter_index_for_ilse, backup_for_ilse, obj_backup_for_ilse, con_backup_for_ilse );
          }
          if( (j == 0) || betterFitness( *obj, *con, obj_ilse_best, con_ilse_best ) || equalFitness( *obj, *con, obj_ilse_best, con_ilse_best ) )
          {
            j_ilse_best   = j;
            obj_ilse_best = *obj;
            con_ilse_best = *con;
          }
          backup_for_ilse[FOSs[gomea_index][i][j]] = populations[gomea_index][index][FOSs[gomea_index][i][j]];
          obj_backup_for_ilse = *obj;
          con_backup_for_ilse = *con;
        }
        for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
        {
          result[FOSs[gomea_index][i][j]]          = backup[FOSs[gomea_index][i][j]];
          backup_for_ilse[FOSs[gomea_index][i][j]] = backup[FOSs[gomea_index][i][j]];
        }
        for( j = 0; j <= j_ilse_best; j++ )
          result[FOSs[gomea_index][i][j]] = populations[gomea_index][index][FOSs[gomea_index][i][j]];
        *obj = obj_ilse_best;
        *con = con_ilse_best;
      }

      if( betterFitness( *obj, *con, obj_backup, con_backup ) || equalFitness( *obj, *con, obj_backup, con_backup ) )
      {
        // TODO: Only copying the FOS component part is bad if the entire solution has changed
        // If an LS procedure manipulates the keys outside of the FOS, things will go wrong!
        if (eval_performs_ls) {
          for( j = 0; j < number_of_parameters; j++ )
            backup[j] = result[j];
        } else {
          for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
            backup[FOSs[gomea_index][i][j]] = result[FOSs[gomea_index][i][j]];
        }
        if( use_ilse )
        {
          for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
            backup_for_ilse[FOSs[gomea_index][i][j]] = result[FOSs[gomea_index][i][j]];
        }

        obj_backup = *obj;
        con_backup = *con;

        solution_has_changed = 1;
      }
      else
      {
        if (eval_performs_ls) {
          for( j = 0; j < number_of_parameters; j++ )
            result[j] = backup[j];
        } else {
          for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
            result[FOSs[gomea_index][i][j]] = backup[FOSs[gomea_index][i][j]];
        }
        *obj = obj_backup;
        *con = con_backup;
      }
    }
  }

  /* Phase 2 (Forced Improvement): optimal mixing with elitist solution */
  if( (!solution_has_changed) || (no_improvement_stretchs[gomea_index] > (10+10*(log(population_sizes[gomea_index])/log(10)))) )
  {
    shuffleFOSSpecificGOMEA( gomea_index );
    if( use_ilse )
      shuffleFOSSubsetsSpecificGOMEA( gomea_index );

    solution_has_changed = 0;
    for( i = 0; i < FOSs_length[gomea_index]; i++ )
    {
      /* Convert elite solution to binary representation and set factor variables. */
      for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
        result[FOSs[gomea_index][i][j]] = elitist_solution[FOSs[gomea_index][i][j]];

      /* Test if the change is for the better */
      donor_parameters_are_the_same = 1;
      for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
      {
        if( backup[FOSs[gomea_index][i][j]] != result[FOSs[gomea_index][i][j]] )
        {
          donor_parameters_are_the_same = 0;
          break;
        }
      }
      if( !donor_parameters_are_the_same )
      {
        if( !use_ilse )
          installedProblemEvaluation( problem_index, result, obj, con, FOSs_number_of_indices[gomea_index][i], FOSs[gomea_index][i], backup, obj_backup, con_backup );
        else
        {
          for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
            result[FOSs[gomea_index][i][j]] = backup[FOSs[gomea_index][i][j]];

          j_ilse_best         = 0;
          obj_ilse_best       = *obj;
          con_ilse_best       = *con;
          obj_backup_for_ilse = *obj;
          con_backup_for_ilse = *con;
          for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
          {
            if( result[FOSs[gomea_index][i][j]] != elitist_solution[FOSs[gomea_index][i][j]] )
            {
              result[FOSs[gomea_index][i][j]] = elitist_solution[FOSs[gomea_index][i][j]];
              parameter_index_for_ilse[0] = FOSs[gomea_index][i][j];
              installedProblemEvaluation( problem_index, result, obj, con, 1, parameter_index_for_ilse, backup_for_ilse, obj_backup_for_ilse, con_backup_for_ilse );
            }
            if( (j == 0) || betterFitness( *obj, *con, obj_ilse_best, con_ilse_best ) || equalFitness( *obj, *con, obj_ilse_best, con_ilse_best ) )
            {
              j_ilse_best   = j;
              obj_ilse_best = *obj;
              con_ilse_best = *con;
            }
            backup_for_ilse[FOSs[gomea_index][i][j]] = elitist_solution[FOSs[gomea_index][i][j]];
            obj_backup_for_ilse = *obj;
            con_backup_for_ilse = *con;
          }
          for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
          {
            result[FOSs[gomea_index][i][j]]          = backup[FOSs[gomea_index][i][j]];
            backup_for_ilse[FOSs[gomea_index][i][j]] = backup[FOSs[gomea_index][i][j]];
          }
          for( j = 0; j <= j_ilse_best; j++ )
            result[FOSs[gomea_index][i][j]] = elitist_solution[FOSs[gomea_index][i][j]];
          *obj = obj_ilse_best;
          *con = con_ilse_best;
        }

        if( betterFitness( *obj, *con, obj_backup, con_backup ) )
        {
          if (eval_performs_ls) {
            for( j = 0; j < number_of_parameters; j++ )
              backup[j] = result[j];
          } else {
            for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
              backup[FOSs[gomea_index][i][j]] = result[FOSs[gomea_index][i][j]];
          }
          if( use_ilse )
          {
            for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
              backup_for_ilse[FOSs[gomea_index][i][j]] = result[FOSs[gomea_index][i][j]];
          }

          obj_backup = *obj;
          con_backup = *con;

          solution_has_changed = 1;
        }
        else
        {
          if (eval_performs_ls) {
            for( j = 0; j < number_of_parameters; j++ )
              result[j] = backup[j];
          } else {
            for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
              result[FOSs[gomea_index][i][j]] = backup[FOSs[gomea_index][i][j]];
          }

          *obj = obj_backup;
          *con = con_backup;
        }
      }
      if( solution_has_changed )
        break;
    }

    if( !solution_has_changed )
    {
      if( betterFitness( elitist_solution_objective_value, elitist_solution_constraint_value, *obj, *con ) )
        solution_has_changed = 1;

      for( i = 0; i < number_of_parameters; i++ )
        result[i] = elitist_solution[i];
      *obj = elitist_solution_objective_value;
      *con = elitist_solution_constraint_value;
    }
  }

  free( backup );
  if( use_ilse )
    free( backup_for_ilse );

  return( result );
}

/**
 * Shuffles the FOS (ordering), but not the contents
 * of the linkage subsets themselves.
 */
void shuffleFOSSpecificGOMEA( int gomea_index )
{
  int i, *order, **FOSs_new, *FOSs_number_of_indices_new;

  FOSs_new                   = (int **) Malloc( FOSs_length[gomea_index]*sizeof( int * ) );
  FOSs_number_of_indices_new = (int *) Malloc( FOSs_length[gomea_index]*sizeof( int ) );
  order                     = randomPermutation( FOSs_length[gomea_index] );
  for( i = 0; i < FOSs_length[gomea_index]; i++ )
  {
    FOSs_new[i]                   = FOSs[gomea_index][order[i]];
    FOSs_number_of_indices_new[i] = FOSs_number_of_indices[gomea_index][order[i]];
  }
  free( FOSs[gomea_index] );
  free( FOSs_number_of_indices[gomea_index] );
  FOSs[gomea_index]                   = FOSs_new;
  FOSs_number_of_indices[gomea_index] = FOSs_number_of_indices_new;

  free( order );
}

/**
 * Shuffles the linkage subsets (ordering) in the FOS, but not the FOS itself.
 */
void shuffleFOSSubsetsSpecificGOMEA( int gomea_index )
{
  int i, j, *order, *FOSs_subset_new;

  for( i = 0; i < FOSs_length[gomea_index]; i++ )
  {
    order = randomPermutation( FOSs_number_of_indices[gomea_index][i] );

    FOSs_subset_new = (int *) Malloc( FOSs_number_of_indices[gomea_index][i]*sizeof( int ) );
    for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
      FOSs_subset_new[j] = FOSs[gomea_index][i][order[j]];
    free( FOSs[gomea_index][i] );
    FOSs[gomea_index][i] = FOSs_subset_new;

    free( order );
  }
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Ezilaitini -=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * Undoes GOMEA initializations.
 */
void ezilaitiniAllGOMEAs()
{
  int i;

  for( i = 0; i < number_of_GOMEAs; i++ )
    ezilaitiniSpecificGOMEA( i );

  free( FOSs_length );
  free( FOSs_number_of_indices );
  free( FOSs );
  free( dependency_matrices );
  free( populations );
  free( objective_values );
  free( constraint_values );
  free( average_objective_values );
  free( average_constraint_values );
  free( terminated );
  free( offsprings );
  free( objective_values_offsprings );
  free( constraint_values_offsprings );
  free( objective_values_best_of_generation );
  free( constraint_values_best_of_generation );
}

void ezilaitiniSpecificGOMEA( int gomea_index )
{
  int i;

  if( FOSs[gomea_index] != NULL )
  {
    for( i = 0; i < FOSs_length[gomea_index]; i++ )
      free( FOSs[gomea_index][i] );
    free( FOSs[gomea_index] );
    free( FOSs_number_of_indices[gomea_index] );
  }

  for( i = 0; i < number_of_parameters; i++ )
    free( dependency_matrices[gomea_index][i] );
  free( dependency_matrices[gomea_index] );

  ezilaitiniSpecificGOMEAMemoryForPopulationAndOffspring( gomea_index );
}

/**
 * Initializes the memory for a single GOMEA thread, for the population only.
 */
void ezilaitiniSpecificGOMEAMemoryForPopulationAndOffspring( int gomea_index )
{
  int i;

  for( i = 0; i < population_sizes[gomea_index]; i++ )
    free( offsprings[gomea_index][i] );

  for( i = 0; i < population_sizes[gomea_index]; i++ )
    free( populations[gomea_index][i] );

  free( populations[gomea_index] );
  free( objective_values[gomea_index] );
  free( constraint_values[gomea_index] );
  free( offsprings[gomea_index] );
  free( objective_values_offsprings[gomea_index] );
  free( constraint_values_offsprings[gomea_index] );
}

void ezilaitiniValueAndSetOfSolutionsToReach()
{
  int i;

  for( i = 0; i < number_of_solutions_in_sostr; i++ )
    free( sostr[i] );
  free( sostr );
}

void ezilaitiniProblem( int index )
{
  // Nothing to clean up: that is the job of the calling environment!
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Time -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
long getMilliSecondsRunning()
{
  return( getMilliSecondsRunningSinceTimeStamp( timestamp_start ) );
}

long getMilliSecondsRunningAfterInit()
{
  return( getMilliSecondsRunningSinceTimeStamp( timestamp_start_after_init ) );
}

long getMilliSecondsRunningSinceTimeStamp( long timestamp )
{
  long timestamp_now, difference;

  timestamp_now = getCurrentTimeStampInMilliSeconds();

  difference = timestamp_now-timestamp;

  return( difference );
}

long getCurrentTimeStampInMilliSeconds()
{
  struct timeval tv;
  long   result;

  gettimeofday( &tv, NULL );
  result = (tv.tv_sec * 1000) + (tv.tv_usec / 1000);

  return( result );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Run -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Initializes the VOSOSTR, the random number generator and the problem and runs the GOMEA.
 */
void run()
{
  timestamp_start = getCurrentTimeStampInMilliSeconds();

  initializeValueAndSetOfSolutionsToReach();

  initializeRandomNumberGenerator();

  // if( print_verbose_overview )
  //   printVerboseOverview();

  timestamp_start_after_init = getCurrentTimeStampInMilliSeconds();

  multiPopGOMEA();

  // writeRunningTime( (char *) "total_running_time.dat" );

  ezilaitiniProblem( problem_index );

  ezilaitiniValueAndSetOfSolutionsToReach();
}

void multiPopGOMEA()
{
  maximum_number_of_GOMEAs                  = 25;
  number_of_subgenerations_per_GOMEA_factor = 4;
  base_population_size                      = 1;

  number_of_GOMEAs               = 0;
  number_of_generations          = 0;
  number_of_evaluations          = 0;
  number_of_bit_flip_evaluations = 0;
  minimum_GOMEA_index            = 0;

  population_sizes                   = (int *) Malloc( maximum_number_of_GOMEAs*sizeof( int ) );
  number_of_subgenerations_per_GOMEA = (int *) Malloc( maximum_number_of_GOMEAs*sizeof( int ) );
  no_improvement_stretchs            = (int *) Malloc( maximum_number_of_GOMEAs*sizeof( int ) );

  while( !checkTermination() )
  {
    if( number_of_GOMEAs < maximum_number_of_GOMEAs )
      initializeNewGOMEA();
    
    if( earlystop )
      break;
    
    if( write_generational_statistics )
      writeGenerationalStatistics();

    if( write_generational_solutions )
      writeGenerationalSolutions( 0 );

    generationalStepAllGOMEAs();

    number_of_generations++;

    if( earlystop )
      break;
  }

  if( write_generational_statistics )
    writeGenerationalStatistics();

  if( write_generational_solutions )
    writeGenerationalSolutions( 1 );

  ezilaitiniAllGOMEAs();

  free( no_improvement_stretchs );
  free( number_of_subgenerations_per_GOMEA );
  free( population_sizes );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Main -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * The main function:
 * - interpret parameters on the command line
 * - run the algorithm with the interpreted parameters
 */
// int main( int argc, char **argv )
// {
//   interpretCommandLine( argc, argv );

//   elitist_solution = (double *) Malloc( number_of_parameters*sizeof( double ) );
//   run();
//   free( elitist_solution );

//   return( 0 );
// }
#ifdef __cplusplus
extern "C" { 
#endif
void optimizeGOMEA(void *juliafunctiona, JuliaFunctionRunner calljuliaa, 
  int dim, int fos, int eva, int mil, 
  int randomkeysuntila, char reencodea, 
  char rescalea,
  char translatea,
  char range_type, double *range_lb, double *range_ub,
  double *objective, double *constraint, double *best_solution,
  char eval_performs_ls_a,
  char dependency_matrix_method_a)
{
  juliafunction = juliafunctiona;
  calljulia = calljuliaa;
  // TODO: remove: problem_index?
  number_of_parameters           = dim;
  FOSs_structure_index           = fos;
  maximum_number_of_evaluations  = eva;
  maximum_number_of_milliseconds = mil;
  reencode                       = reencodea;
  randomkeysuntil                  = randomkeysuntila;
  rescale                        = rescalea;
  eval_performs_ls               = eval_performs_ls_a;
  dependency_matrix_method       = dependency_matrix_method_a;
  translate                      = translatea;

  new_solution_generator = range_type;
  if (range_type == 0)
  {
    rangeslb = (double*) malloc(sizeof(double));
    rangesub = (double*) malloc(sizeof(double));
    rangeslb[0] = range_lb[0];
    rangesub[0] = range_ub[0];
  }
  else
  {
    rangeslb = (double *) Malloc(sizeof(double) * dim);
    rangesub = (double *) Malloc(sizeof(double) * dim);
    for( int j = 0; j < dim; j++ )
    {
      rangeslb[j] = range_lb[j];
      rangesub[j] = range_ub[j];
    }
  }

  write_generational_statistics = 0;
  write_generational_solutions  = 0;
  print_verbose_overview        = 0;
  print_FOSs_contents           = 0;
  use_ilse                      = 0;

  earlystop = 0;
  
  elitist_solution = (double *) Malloc( number_of_parameters*sizeof( double ) );
  run();
  // Copy over the solution to a julia array.
  for(int i = 0; i < number_of_parameters; i++)
    // Note, best solution is a julia thing, so this is a bit of a hack.
    best_solution[sizeof(double)+i] = elitist_solution[i]; 
  free( elitist_solution );
  free( rangeslb );
  free( rangesub );

  // if(earlystop) 
  //   printf("Stopping early\n");
  
  // Update final objective and constraint value.
  *objective = elitist_solution_objective_value;
  *constraint = elitist_solution_constraint_value;
  // elitist solution is updated automatically
}
#ifdef __cplusplus
}
#endif

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/