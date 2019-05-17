/*********************                                                        */
/*! \file Reluplex.h
 ** \verbatim
 ** Top contributors (to current version):
 **   Guy Katz
 ** This file is part of the Reluplex project.
 ** Copyright (c) 2016-2017 by the authors listed in the file AUTHORS
 ** (in the top-level source directory) and their institutional affiliations.
 ** All rights reserved. See the file COPYING in the top-level source
 ** directory for licensing information.\endverbatim
 **/

#ifndef __Reluplex_h__
#define __Reluplex_h__

#include "Debug.h"
#include "File.h"
#include "FloatUtils.h"
#include "GlpkWrapper.h"
#include "IReluplex.h"
#include "Map.h"
#include "Queue.h"
#include "ReluPairs.h"
#include "Set.h"
#include "Tableau.h"
#include "Stack.h"
#include "SmtCore.h"
#include "Stringf.h"
#include "Time.h"
#include "VariableBound.h"
#include <string.h>

static const double ALMOST_BROKEN_RELU_MARGIN = 0.001;
static const double GLPK_IMPRECISION_TOLERANCE = 0.001;
static const double NUMBERICAL_INSTABILITY_CONSTANT = 0.0001;
static const double OOB_EPSILON = 0.001;
static const double MAX_ALLOWED_DEGRADATION = 0.000001;

// How often should the statistics and assignment be printed
static const unsigned PRINT_STATISTICS = 500;
static const unsigned PRINT_ASSIGNMENT = 500;

// How many times GLPK is allowed to fail before tableau restoration
static const unsigned MAX_GLPK_FAILURES_BEFORE_RESOTRATION = 10;

class Reluplex;

static Reluplex *activeReluplex;

// Callbacks from GLPK`
void boundCalculationHook( int n, int m, int *head, int leavingBasic, int enteringNonBasic, double *basicRow );
void iterationCountCallback( int count );
void reportSoiCallback( double soi );
int makeReluAdjustmentsCallback( int n, int m, int nonBasicEncoding, const int *head, const char *flags );

class InvariantViolationError
{
public:
    InvariantViolationError( unsigned violatingStackLevel )
            : _violatingStackLevel( violatingStackLevel )
    {
    }

    unsigned _violatingStackLevel;
};

String milliToString( unsigned long long milliseconds )
{
    unsigned seconds = milliseconds / 1000;
    unsigned minutes = seconds / 60;
    unsigned hours = minutes / 60;

    return Stringf( "%02u:%02u:%02u", hours, minutes - hours * 60, seconds - ( minutes * 60 ) );
}

class Reluplex : public IReluplex
{
public:
    /*** Add by lzs **/
    bool stopFind;  // Add by lzs
    double *_assignment;	//数组，存储所有变量的相应值
    int *candidateNode;

    /*** Add end **/

    enum FinalStatus {
        SAT = 0,
        UNSAT = 1,
        ERROR = 2,
        NOT_DONE = 3,
        SAT_BUT_THEN_UNSAT = 4,
        SAT_BUT_THEN_ERROR = 5,
        SAT_BUT_THEN_NOT_DONE = 6,
    };

    struct GlpkRowEntry
    {
        unsigned _variable;
        double _coefficient;

        GlpkRowEntry( unsigned variable, double coefficient )
                : _variable( variable )
                , _coefficient( coefficient )
        {
        }
    };

    Reluplex( unsigned numVariables, char *finalOutputFile = NULL, String reluplexName = "" )
            :   stopFind(false)
            , _assignment( NULL )
            ,candidateNode(NULL)
            ,_numVariables( numVariables )	//所有变量的总数
            , _reluplexName( reluplexName )	//
            , _finalOutputFile( finalOutputFile )
            , _finalStatus( NOT_DONE )
            , _wasInitialized( false )
            , _tableau( numVariables )
            , _preprocessedTableau( numVariables )
            , _upperBounds( NULL )	// 所有和bounds和assignment有关的变量，都是数组，存储所有变量的上下界和赋值
            , _lowerBounds( NULL )
            , _preprocessedUpperBounds( NULL )
            , _preprocessedLowerBounds( NULL )
            , _preprocessedAssignment( NULL )// 所有和bounds和assignment有关的变量，都是数组，存储所有变量的上下界和赋值
            , _smtCore( this, _numVariables )
            , _useApproximations( true )
            , _findAllPivotCandidates( false )
            , _conflictAnalysisCausedPop( 0 )
            , _logging(false )     // change by lzs
            , _dumpStates( false )
            , _numCallsToProgress( 0 )
            , _numPivots( 0 )
            , _totalPivotTimeMilli( 0 )
            , _totalDegradationCheckingTimeMilli( 0 )
            , _totalRestorationTimeMilli( 0 )
            , _totalPivotCalculationCount( 0 )
            , _totalNumBrokenRelues( 0 )
            , _brokenRelusFixed( 0 )
            , _brokenReluFixByUpdate( 0 )
            , _brokenReluFixByPivot( 0 )
            , _brokenReluFixB( 0 )
            , _brokenReluFixF( 0 )
            , _numEliminatedVars( 0 )
            , _varsWithInfiniteBounds( 0 )
            , _numStackSplits( 0 )
            , _numStackMerges( 0 )
            , _numStackPops( 0 )
            , _numStackVisitedStates( 0 )
            , _currentStackDepth( 0 )
            , _minStackSecondPhase( 0 )
            , _maximalStackDepth( 0 )
            , _boundsTightendByTightenAllBounds( 0 )
            , _almostBrokenReluPairCount( 0 )
            , _almostBrokenReluPairFixedCount( 0 )
            , _numBoundsDerivedThroughGlpk( 0 )
            , _numBoundsDerivedThroughGlpkOnSlacks( 0 )
            , _totalTightenAllBoundsTime( 0 )
            , _eliminateAlmostBrokenRelus( false )
            , _printAssignment( false )
            , _numOutOfBoundFixes( 0 )
            , _numOutOfBoundFixesViaBland( 0 )
            , _useDegradationChecking( false )
            , _numLpSolverInvocations( 0 )
            , _numLpSolverFoundSolution( 0 )
            , _numLpSolverNoSolution( 0 )
            , _numLpSolverFailed( 0 )
            , _numLpSolverIncorrectAssignment( 0 )
            , _totalLpSolverTimeMilli( 0 )
            , _totalLpExtractionTime( 0 )
            , _totalLpPivots( 0 )
            , _maxLpSolverTimeMilli( 0 )
            , _numberOfRestorations( 0 )
            , _maxDegradation( 0.0 )
            , _totalProgressTimeMilli( 0 )
            , _timeTighteningGlpkBoundsMilli( 0 )
            , _currentGlpkWrapper( NULL )
            , _relusDissolvedByGlpkBounds( 0 )
            , _glpkSoi( 0 )
            , _storeGlpkBoundTighteningCalls( 0 )
            , _storeGlpkBoundTighteningCallsOnSlacks( 0 )
            , _storeGlpkBoundTighteningIgnored( 0 )
            , _maxBrokenReluAfterGlpk( 0 )
            , _totalBrokenReluAfterGlpk( 0 )
            , _totalBrokenNonBasicReluAfterGlpk( 0 )
            , _useSlackVariablesForRelus( USE_ROW_SLACK_VARIABLES )
            , _fixRelusInGlpkAssignmentFixes( 0 )
            , _fixRelusInGlpkAssignmentInvoked( 0 )
            , _fixRelusInGlpkAssignmentIgnore( 0 )
            , _maximalGlpkBoundTightening( false )
            , _useConflictAnalysis( true )
            , _temporarilyDontUseSlacks( false )
            , _quit( false )
            , _fullTightenAllBounds( true )
            , _glpkExtractJustBasics( true )
            , _totalTimeEvalutingGlpkRows( 0 )
            , _consecutiveGlpkFailureCount( 0 )
            , alreadySAT(false)
    {

        activeReluplex = this;

        srand( time( 0 ) );

        _upperBounds = new VariableBound[_numVariables];
        _lowerBounds = new VariableBound[_numVariables];
        _preprocessedUpperBounds = new VariableBound[_numVariables];
        _preprocessedLowerBounds = new VariableBound[_numVariables];

        _assignment = new double[_numVariables];
        _preprocessedAssignment = new double[_numVariables];

        for ( unsigned i = 0; i < _numVariables; ++i )
        {
            _assignment[i] = 0.0;
            _preprocessedAssignment[i] = 0.0;
        }

        FloatUtils::printEpsion();
        printf( "Almost-broken nuking marging: %.15lf\n", ALMOST_BROKEN_RELU_MARGIN );
    }

    ~Reluplex()
    {
        if ( _finalOutputFile )
        {
            printFinalStatistics();
        }

        if ( _upperBounds )
        {
            delete[] _upperBounds;
            _upperBounds = NULL;
        }

        if ( _lowerBounds )
        {
            delete[] _lowerBounds;
            _lowerBounds = NULL;
        }

        if ( _preprocessedLowerBounds )
        {
            delete[] _preprocessedLowerBounds;
            _preprocessedLowerBounds = NULL;
        }

        if ( _preprocessedUpperBounds )
        {
            delete[] _preprocessedUpperBounds;
            _preprocessedUpperBounds = NULL;
        }

        if ( _assignment )
        {
            delete[] _assignment;
            _assignment = NULL;
        }

        if ( _preprocessedAssignment )
        {
            delete[] _preprocessedAssignment;
            _preprocessedAssignment = NULL;
        }
    }


    void initialize()
    {
        // 检查所有变量的上下界是否是合理值（上界大于下界）
        // 以及非基变量non-basic是否在上下界之内

        //（此时的non-basic是所有具有真实意义的值，因为此时没有经过pivot，所有人工变量都是basic）
        // 如果non-baisc有越界值，要更新值使非基变量non-basic在上下界内，（即使使basic越界也无所谓
        // 如果被update的值是relu中的一个，要立即考虑是否需要fix，如果需要，立即在update中用updateB或updateF来进行fix

        /*** make all non-baisc in their bounds ***/
        // 将现阶段所有的non-basic都update到上下界之内，没有执行过Pivot
        initialUpdate();

        // 处理无穷大变量，以及non-basic的上下界tighter
        makeAllBoundsFinite();

        _wasInitialized = true;
    }

    FinalStatus solve()
    {
        timeval start = Time::sampleMicro();
        timeval end;

        try
        {
            if ( !_wasInitialized )
                initialize();

            countVarsWithInfiniteBounds();
            if ( !eliminateAuxVariables() )
            {
                _finalStatus = Reluplex::ERROR;
                end = Time::sampleMicro();
                _totalProgressTimeMilli += Time::timePassed( start, end );
                return _finalStatus;
            }

            storePreprocessedMatrix();

            printf( "Initialization steps over.\n" );
            printStatistics();
            dump();
            printf( "Starting the main loop\n" );

            while ( !_quit )
            {
                computeVariableStatus();

                if ( allVarsWithinBounds() && allRelusHold() )
                {
                    dump();
                    showDissolvedMergeReluPairs();
                    printStatistics();
                    _finalStatus = Reluplex::SAT;
                    end = Time::sampleMicro();
                    _totalProgressTimeMilli += Time::timePassed( start, end );
                    return _finalStatus;
                }

                unsigned violatingLevelInStack;
                if ( !progress( violatingLevelInStack ) )
                {
                    if ( _useConflictAnalysis )
                        _smtCore.pop( violatingLevelInStack );
                    else
                        _smtCore.pop();

                    setMinStackSecondPhase( _currentStackDepth );
                }
            }
        }
        catch ( const Error &e )
        {
            end = Time::sampleMicro();
            _totalProgressTimeMilli += Time::timePassed( start, end );

            if ( e.code() == Error::STACK_IS_EMPTY )
            {
                _finalStatus = Reluplex::UNSAT;
                return _finalStatus;
            }
            else
            {
                printf( "Found error: %u\n", e.code() );

                _finalStatus = Reluplex::ERROR;
                return _finalStatus;
            }
        }
        catch ( const InvariantViolationError &e )
        {
            end = Time::sampleMicro();
            _totalProgressTimeMilli += Time::timePassed( start, end );
            _finalStatus = Reluplex::UNSAT;
            return _finalStatus;
        }
        catch ( ... )
        {
            end = Time::sampleMicro();
            _totalProgressTimeMilli += Time::timePassed( start, end );
            _finalStatus = Reluplex::ERROR;
            return _finalStatus;
        }

        // Quit was called
        _finalStatus = Reluplex::NOT_DONE;
        end = Time::sampleMicro();
        _totalProgressTimeMilli += Time::timePassed( start, end );
        return _finalStatus;
    }

    void printBounds(){
//        printf("\n~~~~~printBounds()\n");
        for (unsigned i = 0; i < _numVariables; i++) {
            printf("the var : %u \t", i);
            if ( _lowerBounds[i].finite() ){
                printf("_lowerBounds : %5.2lf \t", _lowerBounds[i].getBound());
            } else{
                printf("infinite \t");
            }

            if ( _upperBounds[i].finite() ){
                printf("_upperBounds : %5.2lf \t", _upperBounds[i].getBound());
            }
            else{
                printf("infinite\t");
            }
            printf("\n");
        }
    }


    FinalStatus solve(double **currentAdversaryE, unsigned &num_AE, unsigned &num_Node, unsigned &num_Expected_AE )
    {
        timeval start = Time::sampleMicro();
        timeval end;

        printBounds();

        try
        {
            /*** first step*****/
            // 处理non-basic变量的越界问题，以及对无穷界限进行缩小
            if ( !_wasInitialized )
                initialize();

            countVarsWithInfiniteBounds();

            /*** second step*****/
            // 消除所有另外引入的辅助变量,即所有basic变量，即使只有一个无法消除，也会使整个项目failed
            if ( !eliminateAuxVariables() )
            {
                _finalStatus = Reluplex::ERROR;

                //// add by lzs
                if (alreadySAT){
                    _finalStatus = Reluplex::SAT_BUT_THEN_ERROR;
                }
                //// add end
                end = Time::sampleMicro();
                _totalProgressTimeMilli += Time::timePassed( start, end );
                return _finalStatus;
            }


            /** third step*****/
            // 保存一些需要的数据，如initialize后的Tableau、全文只在这里有调用
            storePreprocessedMatrix();

            // 打印信息和存入Log
            printf( "Initialization steps over.\n" );
            printf("\n----- printStatistics():after initialize() ----\n");
            printStatistics();
            dump();
            printf( "Starting the main loop\n" );


            /** forth step*****/
            while ( !_quit )
            {
                computeVariableStatus();
                /** fifth step*****/
                // 如果满足这两个条件，则返回SAT
                if ( allVarsWithinBounds() && allRelusHold() )
                {
                    //// add by lzs : 在满足条件的时候进行一次判断，要求Basic里面没有只有本身为-1的情况，否则这样的情况就算满足条件，也一般不满足，所以要去除，回退上一步
                    if (checkIfJustItself()) {
                        _smtCore.pop();
                        setMinStackSecondPhase( _currentStackDepth );   // 模仿原代码
                        continue;   // 还原上一步之后，要跳过后面的代码，继续循环
                    }
                    //// add end

                    log( "\nIt can be solved. current state is: \n" );
                    dump();
                    showDissolvedMergeReluPairs();
                    printf("\n----- printStatistics():find an SAT answer ----\n");
                    printStatistics();
                    _finalStatus = Reluplex::SAT;
                    end = Time::sampleMicro();
                    _totalProgressTimeMilli += Time::timePassed( start, end );

                    /*** add by lzs : adv ***/
                    alreadySAT = true;
                    Map<unsigned, unsigned> indexToVar;
                    /******** change! **********/
                    // cto
//                    unsigned inputLayerSize = 18;   // 输入层有多少个元素
//                    unsigned outputLayerSize = 3;   // 输出层有多少个元素
//                    unsigned baseOutputIndex = 234; // 输出层第一个元素的下标
                    // cto end

                    // leaky_example_new
                    unsigned inputLayerSize = 1;   // 输入层有多少个元素
                    unsigned outputLayerSize = 1;   // 输出层有多少个元素
                    unsigned baseOutputIndex = 5; // 输出层第一个元素的下标
                    // leaky_example_new end

                    // 其实以后只需要改上面的参数就好
                    for (unsigned ipt = 0; ipt < inputLayerSize; ipt++) {
                        indexToVar[ipt] = ipt;
                    }
                    for (unsigned otpt = 0; otpt < outputLayerSize; ++otpt) {
                        indexToVar[inputLayerSize + otpt] = baseOutputIndex + otpt;
                    }
                    /******* change end ******/

                    for (unsigned c = 0; c < num_Node; c++) {
                        currentAdversaryE[num_AE][c] = _assignment[ indexToVar[c] ];
                    }
                    num_AE ++;  //表示当前找到的个数，也表示下一个应该存储的下标值

                    printf("printCurrentAE in solve()'s main loop:\n");
//                    for (unsigned j = 0; j < num_AE; ++j) {   // 输出AE数组里的全部内容
                    for (unsigned j = num_AE - 1; j < num_AE; ++j) {     //简便输出，只输出最近找到的一个
                        printf("This is a adversial example: %u \n", num_AE);
                        for (unsigned k = 0; k < num_Node; ++k) {
                            if( k < inputLayerSize){
                                printf( "Input Variable %u : value = %.10lf \n", k , currentAdversaryE[j][k]);
                            } else{
                                printf( "Output Variable %u : value = %.10lf \n", k , currentAdversaryE[j][k]);
                            }
                        }
                    }

                    if (num_AE < num_Expected_AE){
                        printf("\nWe find a AE, but will run continuely: %u \n",num_AE);
                        _smtCore.pop();
                        setMinStackSecondPhase( _currentStackDepth );   // 模仿原代码
                        printf("After pop(), current assignment is : \n");
                        for (unsigned c = 0; c < num_Node; c++) {
                            if(c < inputLayerSize){  //
                                printf("input[%u]:%.10lf \n",c,_assignment[ indexToVar[c] ]);
                            } else{
                                printf("output[%u]:%.10lf \n",c,_assignment[ indexToVar[c] ]);
                            }
                        }
                        printf("\n");
                        continue;
                    } else{
                        return _finalStatus;
                    }
                    /*** add end ***/
                }
                // violatingLevelInStack默认初始化为0，表示导致violation进行了多少次split,而smtCore需要根据这个值重做多少次decisions
                unsigned violatingLevelInStack;

                // 否则，还有越界的basic变量，或者broken的relu，就调用process继续处理
                // 如果处理成功，表明进行了一次值的更新，返回true,进行下一次循环

                // 如果处理不成功，返回false
                // 可能导致返回false的情况：1、GLPK没有找到解决方案NO_SOLUTION_EXISTS
                // 2、出现异常
                // 3、fix越界变量时没有找到可以Pivot的候选者

                // 此时调用smtScore，根据_useConflictAnalysis，对栈进行pop一次或者多次，pop之后之前的状态直接就还原进入了_reluplex对象
                // 进入下一次循环，继续处理

                //// 每次progress可能是用GLPK进行一次fix越界变量，或者是用update-b-f更新有问题的pair，或者是用smtCore进行split,并记录状态在stack中

                if ( !progress( violatingLevelInStack ) )
                {
                    // 如果需要进行冲突处理，则调用_smtCore, default value is true
                    // violatingLevelInStack的值是刚进入progress时的_currentStackDepth，

                    if ( _useConflictAnalysis )
                        _smtCore.pop( violatingLevelInStack );  // if has conflict Analysis, it will go back several steps
                    else
                        _smtCore.pop();     // if has no , it will just go back one step

                    setMinStackSecondPhase( _currentStackDepth );
                }

                // 如果progress处理成功，返回true,可能是SOLVER_FAILED，也可能是SOLVER_FOUND，
                // 进入下个循环

            }
        }
        catch ( const Error &e )
        {
            end = Time::sampleMicro();
            _totalProgressTimeMilli += Time::timePassed( start, end );

            if ( e.code() == Error::STACK_IS_EMPTY )
            {
                _finalStatus = Reluplex::UNSAT;
                //// add by lzs
                if( alreadySAT ) {
                    _finalStatus = Reluplex::SAT_BUT_THEN_UNSAT;
                }
                //// add end
                return _finalStatus;
            }
            else
            {
                printf( "Found error: %u\n", e.code() );
                _finalStatus = Reluplex::ERROR;
                //// add by lzs
                if( alreadySAT ) {
                    _finalStatus = Reluplex::SAT_BUT_THEN_ERROR;
                }
                //// add end
                return _finalStatus;
            }
        }
        catch ( const InvariantViolationError &e )
        {
            end = Time::sampleMicro();
            _totalProgressTimeMilli += Time::timePassed( start, end );
            _finalStatus = Reluplex::UNSAT;
            //// add by lzs
            if( alreadySAT ){
                _finalStatus = Reluplex::SAT_BUT_THEN_UNSAT;
            }
            //// add end
            return _finalStatus;
        }
        catch ( ... )
        {
            end = Time::sampleMicro();
            _totalProgressTimeMilli += Time::timePassed( start, end );
            _finalStatus = Reluplex::ERROR;
            //// add by lzs
            if( alreadySAT ){
                _finalStatus = Reluplex::SAT_BUT_THEN_ERROR;
            }
            //// add end
            return _finalStatus;
        }

        // Quit was called
        _finalStatus = Reluplex::NOT_DONE;
        //// add by lzs
        if( alreadySAT ){
            _finalStatus = Reluplex::SAT_BUT_THEN_NOT_DONE;
        }
        //// add end
        end = Time::sampleMicro();
        _totalProgressTimeMilli += Time::timePassed( start, end );
        return _finalStatus;
    }


    bool progress( unsigned &violatingLevelInStack )
    {
        /**
         * progress()有两个出口返回true
         * 1. basic变量中没有越界变量，但整体有broken relu Pairs，成功修复后，返回true
         * 2. basic变量中有越界变量，调用Glpk进行求解，得到SOLVER_FAILED的结果。返回true
         * 3. basic变量中有越界变量，调用Glpk进行求解，得到SOLUTION_FOUND结果，且此时没有broken的relu Pair,返回true.
         * 如果有relu pair问题，要先fix,然后才返回true
         */
        log( "Progress starting\n" );
//        printf( "\n~~~~~~~ _dissolvedReluVariables's size: %u \n", _dissolvedReluVariables.size() );

        try
        {
            ++_numCallsToProgress;

            // The default
            violatingLevelInStack = _currentStackDepth;  //堆栈深度，进行了多少次split

            // _useDegradationChecking不知道什么鬼，默认为false,先不管
            if ( _useDegradationChecking && ( _numCallsToProgress % 50 == 0 ) )
            {
                double currentMaxDegradation = checkDegradation();
                if ( currentMaxDegradation > MAX_ALLOWED_DEGRADATION )
                {
                    restoreTableauFromBackup();
                    return true;
                }
            }

            // 每调用500次进行一次数据打印
            if ( _numCallsToProgress % PRINT_STATISTICS == 0 )
                printStatistics();

            // 每500次打印赋值
            if ( _printAssignment && _numCallsToProgress % PRINT_ASSIGNMENT == 0 )
                printAssignment();

            dump();

            /********以下开始papaer中所写**********/

            List<unsigned> outOfBoundVariables;
            findOutOfBounds( outOfBoundVariables ); //在basic中找到越界变量，存入List中

            // If we have out-of-bounds variables, we deal with them first
            if ( !outOfBoundVariables.empty() )
            {
                log( "Progress: have OOB vars\n" );
//                printf("\n~~~~~glpk begin\n");

                GlpkWrapper::GlpkAnswer answer = fixOutOfBounds();

                // _consecutiveGlpkFailureCount用于记录GLPK Failures的次数=10，如果大于，则报错
                if ( _consecutiveGlpkFailureCount > MAX_GLPK_FAILURES_BEFORE_RESOTRATION )
                {
                    printf( "Error: %u Consecutive GLPK failures\n", MAX_GLPK_FAILURES_BEFORE_RESOTRATION );
                    throw Error( Error::CONSECUTIVE_GLPK_FAILURES );
                }

                // 如果返回 NO_SOLUTION_EXISTS，表示没有找到解决方案
                if ( answer == GlpkWrapper::NO_SOLUTION_EXISTS )
                    return false;

                // 如果返回 SOLVER_FAILED，表示本次解决失败，回溯到之前的步骤再解决
                if ( answer == GlpkWrapper::SOLVER_FAILED )
                {
                    // In this case, we restored from the original tableau; nothing left to do here.
                    return true;
                }

                //// 如果返回值非以上两种，则只可能是 SOLUTION_FOUND，表示找到了解决方案
                // 此时如果没有broken的ReluPair，则直接返回true
                if ( allRelusHold() ){
//                    printf("\n~~~~~after glpk and all the relu hold\n");
                    return true;
                }

//                 Glpk solved, but we have a relu problem. See if any relu pairs can be eliminated.
//                if ( learnedGlpkBounds() )
//                {
//                    timeval boundStart = Time::sampleMicro();
//
//                    unsigned numDissolvedRelusBefore = countDissolvedReluPairs();
//                    try
//                    {
//                        printf("\n~~~~~~try performGlpkBoundTightening\n");
//                        performGlpkBoundTightening();
//                        tightenAllBounds();
//                    }
//                    catch ( ... )
//                    {
//                        timeval boundEnd = Time::sampleMicro();
//                        _timeTighteningGlpkBoundsMilli += Time::timePassed( boundStart, boundEnd );
//                        throw;
//                    }
//
//                    unsigned numDissolvedRelusAfter = countDissolvedReluPairs();
//
//                    timeval boundEnd = Time::sampleMicro();
//                    _timeTighteningGlpkBoundsMilli += Time::timePassed( boundStart, boundEnd );
//
//                    if ( numDissolvedRelusAfter > numDissolvedRelusBefore )
//                    {
//                        _relusDissolvedByGlpkBounds += ( numDissolvedRelusAfter - numDissolvedRelusBefore );
//                    }
//                }
//                printf("\n~~~~~~leave glpk process and begin next progress\n");
                return true;
            }

            //// 没有越界变量，就要查看是否有broken的pair
            // Reset the GLPK failure measures
            _consecutiveGlpkFailureCount = 0;
            _previousGlpkAnswer = GlpkWrapper::SOLUTION_FOUND;

            log( "No OOB variables to fix, looking at broken relus\n" );

            // If we got here, either there are no OOB variables, or they were fixed without changing the tableau
            // and we still have broken relus. Split on one of them.
            List<unsigned> brokenRelus;
            findBrokenRelues( brokenRelus );
            _totalNumBrokenRelues += brokenRelus.size();

            unsigned brokenReluVar = *brokenRelus.begin();
            // 找到broken的ReluPair中的forward变量
            unsigned f = _reluPairs.isF( brokenReluVar ) ? brokenReluVar : _reluPairs.toPartner( brokenReluVar );

            //// notifyBrokenRelu()如果返回true,表示update-f或update-b已经超过阈值，在notifyBrokenRelu中进行了split，并将状态存入栈中
            if ( _smtCore.notifyBrokenRelu( f ) ){
//                return true; // Splitting/Merging is a form of progress,进行过split，直接返回true，以便进入下一次progress循环
                // return false表示在split负数的过程中找不到合适的candidate，已有的系数和leakyRatio不一致，需要回退这一次split


                // return true表示split成功，可以进行后续的正常求解
            }

            //// 否则，若是notifyBrokenRelu返回false，表示此时还不需要进行split，可以继续用update-f-b来进行修复pair
            return fixBrokenRelu( f );// fix完成一次，返回true
        }

        catch ( const InvariantViolationError &e )
        {
            log( "\n\n*** Upper/lower invariant violated! Failure ***\n\n" );

            if ( e._violatingStackLevel != _currentStackDepth )
            {
                violatingLevelInStack = e._violatingStackLevel;
            }

            return false;
        }
    }
    bool progress_temp( unsigned &violatingLevelInStack )
    {
        /**
         * progress()有两个出口返回true
         * 1. basic变量中没有越界变量，但整体有broken relu Pairs，成功修复后，返回true
         * 2. basic变量中有越界变量，调用Glpk进行求解，得到SOLVER_FAILED的结果。返回true
         * 3. basic变量中有越界变量，调用Glpk进行求解，得到SOLUTION_FOUND结果，且此时没有broken的relu Pair,返回true.
         * 如果有relu pair问题，要先fix,然后才返回true
         */

//        copy_reluplex( myCopyReluplex );


        log( "\nProgress starting\n" );

        try
        {
            ++_numCallsToProgress;

            // The default
            violatingLevelInStack = _currentStackDepth;	//堆栈深度，进行了多少次split

            // _useDegradationChecking不知道什么鬼，默认为false,先不管
            if ( _useDegradationChecking && ( _numCallsToProgress % 50 == 0 ) )
            {
                double currentMaxDegradation = checkDegradation();
                if ( currentMaxDegradation > MAX_ALLOWED_DEGRADATION )
                {
                    restoreTableauFromBackup();
                    return true;
                }
            }

            // 每调用500次进行一次数据打印
            if ( _numCallsToProgress % PRINT_STATISTICS == 0 ){
                printf("\n----- printStatistics():500 ----\n");
                printStatistics();
            }
            // 每500次打印赋值
            if ( _printAssignment && _numCallsToProgress % PRINT_ASSIGNMENT == 0 )
                printAssignment();

            dump();

            /********以下开始papaer中所写**********/

            List<unsigned> outOfBoundVariables;
            findOutOfBounds( outOfBoundVariables );	//在basic中找到越界变量，存入List中

            // If we have out-of-bounds variables, we deal with them first
            if ( !outOfBoundVariables.empty() )
            {
                log( "Progress: have OOB vars\n" );

                GlpkWrapper::GlpkAnswer answer = fixOutOfBounds();

                // _consecutiveGlpkFailureCount用于记录GLPK Failures的次数=10，如果大于，则报错
                if ( _consecutiveGlpkFailureCount > MAX_GLPK_FAILURES_BEFORE_RESOTRATION )
                {
                    printf( "Error: %u Consecutive GLPK failures\n", MAX_GLPK_FAILURES_BEFORE_RESOTRATION );
                    throw Error( Error::CONSECUTIVE_GLPK_FAILURES );
                }

                // 如果返回 NO_SOLUTION_EXISTS，表示没有找到解决方案
                if ( answer == GlpkWrapper::NO_SOLUTION_EXISTS )
                    return false;

                // 如果返回 SOLVER_FAILED，表示本次解决失败，回溯到之前的步骤再解决
                if ( answer == GlpkWrapper::SOLVER_FAILED )
                {
                    // In this case, we restored from the original tableau; nothing left to do here.
                    return true;
                }

                //// 如果返回值非以上两种，则只可能是 SOLUTION_FOUND，表示找到了解决方案
                // 此时如果没有broken的ReluPair，则直接返回true
                if ( allRelusHold() )
                    return true;

                // Glpk solved, but we have a relu problem. See if any relu pairs can be eliminated.

                if ( learnedGlpkBounds() )
                {
                    timeval boundStart = Time::sampleMicro();

                    unsigned numDissolvedRelusBefore = countDissolvedReluPairs();
                    try
                    {
                        performGlpkBoundTightening();
                        tightenAllBounds();
                    }
                    catch ( ... )
                    {
                        timeval boundEnd = Time::sampleMicro();
                        _timeTighteningGlpkBoundsMilli += Time::timePassed( boundStart, boundEnd );
                        throw;
                    }

                    unsigned numDissolvedRelusAfter = countDissolvedReluPairs();

                    timeval boundEnd = Time::sampleMicro();
                    _timeTighteningGlpkBoundsMilli += Time::timePassed( boundStart, boundEnd );

                    if ( numDissolvedRelusAfter > numDissolvedRelusBefore )
                    {
                        _relusDissolvedByGlpkBounds += ( numDissolvedRelusAfter - numDissolvedRelusBefore );
                    }
                }

                return true;
            }

            //// 没有越界变量，就要查看是否有broken的pair

            // Reset the GLPK failure measures
            _consecutiveGlpkFailureCount = 0;
            _previousGlpkAnswer = GlpkWrapper::SOLUTION_FOUND;

            log( "No OOB variables to fix, looking at broken relus\n" );

            // If we got here, either there are no OOB variables, or they were fixed without changing the tableau
            // and we still have broken relus. Split on one of them.
            List<unsigned> brokenRelus;
            findBrokenRelues( brokenRelus );    // 找到broken的ReLU Pair，存储在broeknRelus列表中
            _totalNumBrokenRelues += brokenRelus.size();

            unsigned brokenReluVar = *brokenRelus.begin();
            // 找到broken的ReluPair中的forward变量
            unsigned f = _reluPairs.isF( brokenReluVar ) ? brokenReluVar : _reluPairs.toPartner( brokenReluVar );

            //// notifyBrokenRelu()如果返回true,表示update-f或update-b已经超过阈值，在notifyBrokenRelu中进行了split，并将状态存入栈中
            if ( _smtCore.notifyBrokenRelu( f ) )
                return true; // Splitting/Merging is a form of progress，已进行过split，直接返回true，以便进入下一次progress循环

            //// 否则，若是notifyBrokenRelu返回false，表示此时还不需要进行split，可以继续用update-f-b来进行修复pair
            return fixBrokenRelu( f );  // fix完成一次，返回true
        }

        catch ( const InvariantViolationError &e )
        {
            log( "\n\n*** Upper/lower invariant violated! Failure ***\n\n" );

            if ( e._violatingStackLevel != _currentStackDepth )
            {
                violatingLevelInStack = e._violatingStackLevel;
            }

            return false;
        }
    }

    void toggleDegradationChecking( bool value )
    {
        _useDegradationChecking = value;
    }

    void toggleFullTightenAllBounds( bool value )
    {
        _fullTightenAllBounds = value;
    }

    void toggleGlpkExtractJustBasics( bool value )
    {
        _glpkExtractJustBasics = value;
    }

    void togglePrintAssignment( bool value )
    {
        _printAssignment = value;
    }

    void toggleAlmostBrokenReluEliminiation( bool value )
    {
        if ( value )
            printf( "almost-broken relu elimination has been turned on!\n" );

        _eliminateAlmostBrokenRelus = value;
    }

    bool allVarsWithinBounds( bool print = false ) const
    {
        // Only basic variables can be out-of-bounds

        // 在运算过程中所有的non-basic都是在界限范围内的，所以现在只需要检查basic，如果都在界限范围内，那么就意味着所有变量都在范围内了

        for ( auto i : _basicVariables )
            if ( outOfBounds( i ) )
            {
                if ( print )
                {
                    printf( "Variable %u out of bounds: value = %.10lf, range = [%.10lf, %.10lf]\n",
                            i,
                            _assignment[i],
                            _lowerBounds[i].getBound(),
                            _upperBounds[i].getBound() );
                }

                return false;
            }


        return true;
    }

    bool allRelusHold() const
    {
        const Set<ReluPairs::ReluPair> &pairs( _reluPairs.getPairs() );
        for ( const auto &pair : pairs )
        {
            unsigned b = pair.getB();
            unsigned f = pair.getF();

            if ( ( !_dissolvedReluVariables.exists( f ) ) && reluPairIsBroken( b, f ) )
                return false;
        }

        return true;
    }

    bool reluPairIsBroken( unsigned b, unsigned f ) const
    {
        double bVal = _assignment[b];
        double fVal = _assignment[f];

//        // 当forward变量是0,但backward变量不是0时，返回true表示broken
//        if (FloatUtils::isZero(fVal) && (!FloatUtils::isZero(bVal))) {
//            return true;
//        }
//        // 当forward变量是正数，但与backward变量不相等时，返回true表示broken
//        else if (FloatUtils::isPositive(fVal) && (FloatUtils::areDisequal(fVal, bVal))) {
//            return true;
//        }
//        // 当forward变量是负数，但与backward变量*leakyRatio不相等时，返回true表示broken
//        else if ( !FloatUtils::isPositive(fVal) && (FloatUtils::areDisequal(fVal, leakyRatio * bVal))) {
//            return true;
//        }

        return
                ( FloatUtils::isZero(fVal) && (!FloatUtils::isZero(bVal)) ) ||
                ( FloatUtils::isPositive(fVal) && (FloatUtils::areDisequal(fVal, bVal)) ) ||
                (!FloatUtils::isPositive(fVal) && (FloatUtils::areDisequal(fVal, getLeakyValue() * bVal)));
    }
    bool reluPairIsBroken_temp( unsigned b, unsigned f ) const
    {
        double bVal = _assignment[b];
        double fVal = _assignment[f];

        // 当forward变量是0,但backward变量是正数时，返回true表示broken
        // 当forward变量是正数，但与backward变量不相等时，返回true表示broken

        return
                ( FloatUtils::isZero( fVal ) && FloatUtils::isPositive( bVal ) ) ||
                ( FloatUtils::isPositive( fVal ) && ( FloatUtils::areDisequal( fVal, bVal ) ) );
    }

    unsigned countDissolvedReluPairs() const
    {
        return _dissolvedReluVariables.size();
    }


    unsigned countSplits() const
    {
        unsigned splits = 0;
        for ( const auto &it : _dissolvedReluVariables )
        {
            if ( it.second == TYPE_SPLIT )
                ++splits;
        }

        return splits;
    }

    unsigned countMerges() const
    {
        unsigned merges = 0;

        for ( const auto &it : _dissolvedReluVariables )
        {
            if ( it.second == TYPE_MERGE )
                ++merges;
        }

        return merges;
    }

    unsigned countReluPairsAlmostBroken()
    {
        unsigned count = 0;

        const Set<ReluPairs::ReluPair> &pairs( _reluPairs.getPairs() );
        for ( const auto &pair : pairs )
        {
            unsigned b = pair.getB();
            unsigned f = pair.getF();

            if ( reluPairAlmostBroken( b, f ) )
                ++count;
        }

        return count;
    }

    VariableStatus getVarStatus( unsigned variable ) const
    {
        return _varToStatus.at( variable );
    }

    bool reluPairAlmostBroken( unsigned b, unsigned f ) const
    {
        if ( _dissolvedReluVariables.exists( f ) )
            return false;

        if ( _upperBounds[f].finite() )
        {
            double ub = _upperBounds[f].getBound();
            if ( ( !FloatUtils::isZero( ub ) ) && FloatUtils::lte( ub, ALMOST_BROKEN_RELU_MARGIN ) )
                return true;
        }

        if ( _lowerBounds[b].finite() )
        {
            double lb = _lowerBounds[b].getBound();
            if ( FloatUtils::isNegative( lb ) && FloatUtils::gte( lb, -ALMOST_BROKEN_RELU_MARGIN ) )
                return true;
        }

        return false;
    }

    void findOutOfBounds( List<unsigned> &result ) const
    {
        // Only basic variables can be out-of-bounds
        for ( auto i : _basicVariables )
            if ( outOfBounds( i ) )
                result.append( i );
    }

    void countBrokenReluPairs( unsigned &brokenReluPairs, unsigned &brokenNonBasicReluPairs ) const
    {
        brokenReluPairs = 0;
        brokenNonBasicReluPairs = 0;

        const Set<ReluPairs::ReluPair> &pairs( _reluPairs.getPairs() );
        for ( const auto &pair : pairs )
        {
            unsigned b = pair.getB();
            unsigned f = pair.getF();

            if ( ( !_dissolvedReluVariables.exists( f ) ) && reluPairIsBroken( b, f ) )
            {
                ++brokenReluPairs;
                if ( !_basicVariables.exists( b ) && !_basicVariables.exists( f ) )
                    ++brokenNonBasicReluPairs;
            }
        }
    }

    void findBrokenRelues( List<unsigned> &result ) const
    {
        const Set<ReluPairs::ReluPair> &pairs( _reluPairs.getPairs() );
        for ( const auto &pair : pairs )
        {
            unsigned b = pair.getB();
            unsigned f = pair.getF();

            if ( ( !_dissolvedReluVariables.exists( f ) ) && reluPairIsBroken( b, f ) )
            {
//                printf( "\n~~~~~~~~ find broken relu pairs, b: %u, f: %u\n", b, f );
                result.append( b );
                result.append( f );
            }
        }
    }

    bool partOfBrokenRelu( unsigned variable ) const
    {
        if ( !_reluPairs.isRelu( variable ) )
            return false;

        unsigned partner = _reluPairs.toPartner( variable );
        unsigned b, f;

        if ( _reluPairs.isF( variable ) )
        {
            f = variable;
            b = partner;
        }
        else
        {
            b = variable;
            f = partner;
        }

        return reluPairIsBroken( b, f );
    }

    void countVarsWithInfiniteBounds()
    {
        _varsWithInfiniteBounds = 0;
        for ( unsigned i = 0; i < _numVariables; ++i )
            if ( !_upperBounds[i].finite() || !_lowerBounds[i].finite() )
                ++_varsWithInfiniteBounds;
    }

    void computeVariableStatus()
    {
        for ( unsigned i = 0; i < _numVariables; ++i )
            computeVariableStatus( i );
    }

    void computeVariableStatus( unsigned i )
    {
        double value = _assignment[i];

        if ( _upperBounds[i].finite() && _lowerBounds[i].finite() )
        {
            // Both bounds finite
            double ub = _upperBounds[i].getBound();
            double lb = _lowerBounds[i].getBound();

            if ( FloatUtils::gt( value, ub, OOB_EPSILON ) )
                _varToStatus[i] = VariableStatus::ABOVE_UB;
            else if ( FloatUtils::areEqual( value, ub, OOB_EPSILON ) )
                _varToStatus[i] = FloatUtils::areEqual( lb, ub ) ? VariableStatus::FIXED : VariableStatus::AT_UB;
            else if ( FloatUtils::gt( value, lb, OOB_EPSILON ) )
                _varToStatus[i] = VariableStatus::BETWEEN;
            else if ( FloatUtils::areEqual( value, lb, OOB_EPSILON ) )
                _varToStatus[i] = VariableStatus::AT_LB;
            else
                _varToStatus[i] = VariableStatus::BELOW_LB;
        }

        else if ( !_upperBounds[i].finite() && _lowerBounds[i].finite() )
        {
            // Only lower bound is finite
            double lb = _lowerBounds[i].getBound();

            if ( FloatUtils::gt( value, lb, OOB_EPSILON ) )
                _varToStatus[i] = VariableStatus::BETWEEN;
            else if ( FloatUtils::areEqual( value, lb, OOB_EPSILON ) )
                _varToStatus[i] = VariableStatus::AT_LB;
            else
                _varToStatus[i] = VariableStatus::BELOW_LB;
        }

        else if ( _upperBounds[i].finite() && !_lowerBounds[i].finite() )
        {
            // Only upper bound is finite
            double ub = _upperBounds[i].getBound();

            if ( FloatUtils::gt( value, ub, OOB_EPSILON ) )
                _varToStatus[i] = VariableStatus::ABOVE_UB;
            else if ( FloatUtils::areEqual( value, ub, OOB_EPSILON ) )
                _varToStatus[i] = VariableStatus::AT_UB;
            else
                _varToStatus[i] = VariableStatus::BETWEEN;
        }

        else
        {
            // Both bounds infinite
            _varToStatus[i] = VariableStatus::BETWEEN;
        }
    }

    void printStatistics()
    {
        countVarsWithInfiniteBounds();

        printf( "\n" );
        printf( "%s Statistics update:\n", Time::now().ascii() );
        printf( "\tCalls to 'progress': %u. Total time: %llu milli. Average: %llu milli\n",
                _numCallsToProgress, _totalProgressTimeMilli,
                _numCallsToProgress > 0 ? _totalProgressTimeMilli / _numCallsToProgress : 0 );
        printf( "\tPivot operations: %u. ", _numPivots );
        printf( "Total pivot time: %llu milli.\n", _totalPivotTimeMilli );
        printf( "\tAverage pivot time: %llu milli\n", _numPivots > 0 ? _totalPivotTimeMilli / _numPivots : 0 );
        printf( "\tAverage time per calcuation in pivot: %.5f milli\n",
                _totalPivotCalculationCount > 0 ? ((double)_totalPivotTimeMilli) / _totalPivotCalculationCount : 0 );
        printf( "\tAverage number of calculations in pivot: %llu\n",
                _numPivots > 0 ? _totalPivotCalculationCount / _numPivots : 0 );
        printf( "\tAverage number of broken relues per 'progress': %llu\n",
                _numCallsToProgress > 0 ? _totalNumBrokenRelues / _numCallsToProgress : 0 );
        printf( "\tBroken Relus Fixed: %u (Fs: %u, Bs: %u, fix-by-pivot: %u, fix-by-update: %u)\n",
                _brokenRelusFixed, _brokenReluFixF, _brokenReluFixB, _brokenReluFixByPivot, _brokenReluFixByUpdate );
        printf( "\tRelu-to-OOB step ratio: %u / %u = %lf%%. Avg oob steps per relu: %.02lf.\n",
                _brokenRelusFixed, _numOutOfBoundFixes,
                _numOutOfBoundFixes > 0 ? ((double)_brokenRelusFixed) / _numOutOfBoundFixes : 0,
                _brokenRelusFixed > 0 ? ((double)_numOutOfBoundFixes) / _brokenRelusFixed : 0
        );
        printf( "\tAlmost broken relus encountered: %u. Nuked: %u\n",
                _almostBrokenReluPairCount, _almostBrokenReluPairFixedCount );

        printf( "\tTime in TightenAllBounds: %llu milli. Bounds tightened: %llu\n",
                _totalTightenAllBoundsTime, _boundsTightendByTightenAllBounds );

        printf( "\tRelu pairs dissolved: %u. Num splits: %u. Num merges: %u (remaining: %u / %u)\n",
                _dissolvedReluVariables.size(),
                countSplits(), countMerges(),
                _reluPairs.size() - _dissolvedReluVariables.size(),
                _reluPairs.size());

        printf( "\tNum LP solver invocations: %u. Found solution: %u. No Solution: %u. Failed: %u. "
                "Incorrect assignments: %u.\n",
                _numLpSolverInvocations, _numLpSolverFoundSolution, _numLpSolverNoSolution,
                _numLpSolverFailed, _numLpSolverIncorrectAssignment );
        printf( "\t\tTotal time in LP solver: %llu milli. Max: %u milli. Avg per invocation: %llu milli\n",
                _totalLpSolverTimeMilli,
                _maxLpSolverTimeMilli,
                _numLpSolverInvocations > 0 ? _totalLpSolverTimeMilli / _numLpSolverInvocations : 0 );
        printf( "\t\tNumber of pivots in LP solver: %u. Average time per LP pivot operation: %llu milli\n",
                _totalLpPivots, _totalLpPivots > 0 ? _totalLpSolverTimeMilli / _totalLpPivots : 0 );
        printf( "\t\tTotal time extracting tableaus after LP solved: %llu milli. Average: %llu milli.\n",
                _totalLpExtractionTime, _numLpSolverFoundSolution > 0 ?
                                        _totalLpExtractionTime / _numLpSolverFoundSolution : 0 );
        printf( "\t\tTotal time evaulating GLPK rows: %llu\n", _totalTimeEvalutingGlpkRows );
        printf( "\t\tGlpk bound reports: %llu. On slacks: %llu (= %.lf%%). "
                "Ignored due to small coefficients: %llu. Used: %.2lf%%\n",
                _storeGlpkBoundTighteningCalls,
                _storeGlpkBoundTighteningCallsOnSlacks,
                percents( _storeGlpkBoundTighteningCallsOnSlacks, _storeGlpkBoundTighteningCalls ),
                _storeGlpkBoundTighteningIgnored,
                percents( _storeGlpkBoundTighteningCalls - _storeGlpkBoundTighteningIgnored,
                          _storeGlpkBoundTighteningCalls ) );

        printf( "\t\tNumber of GLPK-derived bounds: %u. On slacks: %u (= %.2lf%%). "
                "Time: %llu milli. Relus consequently dissolved: %u\n",
                _numBoundsDerivedThroughGlpk,
                _numBoundsDerivedThroughGlpkOnSlacks,
                percents( _numBoundsDerivedThroughGlpkOnSlacks, _numBoundsDerivedThroughGlpk ),
                _timeTighteningGlpkBoundsMilli, _relusDissolvedByGlpkBounds );

        printf( "\t\tFix-relu-invariant hook invocations: %llu. Actual repairs: %llu (= %.lf%%). "
                "Ignore to prevent cycles: %llu\n",
                _fixRelusInGlpkAssignmentInvoked, _fixRelusInGlpkAssignmentFixes,
                percents( _fixRelusInGlpkAssignmentFixes, _fixRelusInGlpkAssignmentInvoked ),
                _fixRelusInGlpkAssignmentIgnore );

        printf( "\tAverage number of broken relu pairs after glpk invocation: %lf. Max: %u. "
                "Broken and non-basic pairs: %u\n",
                _numLpSolverFoundSolution > 0 ? ((double)_totalBrokenReluAfterGlpk / _numLpSolverFoundSolution) : 0,
                _maxBrokenReluAfterGlpk,
                _totalBrokenNonBasicReluAfterGlpk );
        printf( "\tVars with infinite bounds: %u / %u\n", _varsWithInfiniteBounds, _numVariables );
        printf( "\tEliminated vars: %u\n", _numEliminatedVars );
        printf( "\tStack: Current depth is: %u (maximal = %u, min second phase = %u).\n"
                "\t       So far: %u splits, %u merges, %u pops. Total visited states: %u\n",
                _currentStackDepth, _maximalStackDepth, _minStackSecondPhase,
                _numStackSplits, _numStackMerges, _numStackPops, _numStackVisitedStates );
        printf( "\t\tPops caused by conflict analysis: %u\n", _conflictAnalysisCausedPop );
        printf( "\t\tTotal time in smtCore: %llu milli\n", _smtCore.getSmtCoreTime() );
        printf( "\tCurrent degradation: %.10lf. Time spent checking: %llu milli. Max measured: %.10lf.\n",
                checkDegradation(), _totalDegradationCheckingTimeMilli, _maxDegradation );
        printf( "\tNumber of restorations: %u. Total time: %llu milli. Average: %lf\n",
                _numberOfRestorations,
                _totalRestorationTimeMilli,
                percents( _totalRestorationTimeMilli, _numberOfRestorations )
        );

        unsigned long long totalUnaccountedFor =
                _totalProgressTimeMilli -
                _totalLpSolverTimeMilli -
                _totalLpExtractionTime -
                _timeTighteningGlpkBoundsMilli -
                _smtCore.getSmtCoreTime() -
                _totalRestorationTimeMilli;

        printf( "\n\n\tSummary: Total: %llu milli\n"
                "\t\t1. GLPK: %llu milli (%.lf%%) \n"
                "\t\t2. Extraction + Postprocessing: %llu milli (%.lf%%)\n"
                "\t\t3. Tightening bounds: %llu milli (%.lf%%)\n"
                "\t\t4. Stack operations: %llu milli (%.lf%%)\n"
                "\t\t5. Tableau restoration operations: %llu milli (%.lf%%)\n"
                "\t\t6. Unaccounted for: %llu milli (%.lf%%)\n",
                _totalProgressTimeMilli,
                _totalLpSolverTimeMilli, percents( _totalLpSolverTimeMilli, _totalProgressTimeMilli ),
                _totalLpExtractionTime, percents( _totalLpExtractionTime, _totalProgressTimeMilli ),
                _timeTighteningGlpkBoundsMilli, percents( _timeTighteningGlpkBoundsMilli, _totalProgressTimeMilli ),
                _smtCore.getSmtCoreTime(), percents( _smtCore.getSmtCoreTime(), _totalProgressTimeMilli ),
                _totalRestorationTimeMilli, percents( _totalRestorationTimeMilli, _totalProgressTimeMilli ),
                totalUnaccountedFor, percents( totalUnaccountedFor, _totalProgressTimeMilli )
        );

        printf( "\n" );
    }

    void printFinalStatistics()
    {
        try
        {
            File outputFile( _finalOutputFile );
            outputFile.open( IFile::MODE_WRITE_APPEND );

            outputFile.write( _reluplexName + ", " );

            String status;
            switch ( _finalStatus )
            {
                case SAT:
                    status = "SAT";
                    break;

                case UNSAT:
                    status = "UNSAT";
                    break;

                case ERROR:
                    status = "ERROR";
                    break;

                case NOT_DONE:
                    status = "NOT_DONE:TIMEOUT";
                    break;

                case SAT_BUT_THEN_UNSAT:
                    status = "SAT_BUT_THEN_UNSAT";
                    break;

                case SAT_BUT_THEN_ERROR:
                    status = "SAT_BUT_THEN_ERROR";
                    break;

                case SAT_BUT_THEN_NOT_DONE:
                    status = "SAT_BUT_THEN_NOT_DONE";
                    break;
            }

            outputFile.write( status + ", " );
            outputFile.write( Stringf( "%llu, %s, ",
                                       _totalProgressTimeMilli, milliToString( _totalProgressTimeMilli ).ascii() ) );
            outputFile.write( Stringf( "%lu, ", _maximalStackDepth ) );
            outputFile.write( Stringf( "%lu\n", _numStackVisitedStates ) );
        }
        catch( ... )
        {
            printf( "Final staitstics printing threw an error!\n" );
        }
    }

    static double percents( double numerator, double denominator )
    {
        if ( FloatUtils::isZero( denominator ) )
            return 0;

        return ( numerator / denominator ) * 100;
    }

    void computeSlackBounds()
    {
        _slackToLowerBound.clear();
        _slackToUpperBound.clear();
        _activeSlackRowVars.clear();
        _activeSlackColVars.clear();

        const Set<ReluPairs::ReluPair> &pairs( _reluPairs.getPairs() );
        for ( const auto &pair : pairs )
        {
            unsigned b = pair.getB();
            unsigned f = pair.getF();

            if ( !_dissolvedReluVariables.exists( f ) )
            {
                unsigned slackRowVar = _fToSlackRowVar[f];
                _activeSlackRowVars.insert( slackRowVar );

                if ( _useSlackVariablesForRelus == USE_ROW_SLACK_VARIABLES )
                {
                    _slackToLowerBound[slackRowVar].setBound( 0.0 );
                    _slackToUpperBound[slackRowVar].setBound( _upperBounds[f].getBound() -
                                                              _lowerBounds[b].getBound() );

                    // The lower bounds exists since level 0.
                    // The upper bounds depend on the current bounds of the variables
                    _slackToLowerBound[slackRowVar].setLevel( 0 );
                    _slackToUpperBound[slackRowVar].setLevel
                            ( std::max( _upperBounds[f].getLevel(), _lowerBounds[b].getLevel() ) );

                }
                else
                {
                    // Use also cols
                    unsigned slackColVar = _fToSlackColVar[f];
                    _activeSlackColVars.insert( slackColVar );

                    _slackToLowerBound[slackColVar].setBound( 0.0 );
                    _slackToUpperBound[slackColVar].setBound( _upperBounds[f].getBound() -
                                                              _lowerBounds[b].getBound() );
                    _slackToUpperBound[slackRowVar].setLevel
                            ( std::max( _upperBounds[f].getLevel(), _lowerBounds[b].getLevel() ) );
                }
            }
        }
    }

    GlpkWrapper::GlpkAnswer fixOutOfBounds()
    {
        ++_numOutOfBoundFixes;
        ++_numLpSolverInvocations;

        timeval lpStart = Time::sampleMicro();
        GlpkWrapper glpkWrapper;
        _currentGlpkWrapper = &glpkWrapper;
        _glpkStoredLowerBounds.clear();
        _glpkStoredUpperBounds.clear();
        _activeSlackRowVars.clear();
        _activeSlackColVars.clear();

        if ( _useSlackVariablesForRelus != DONT_USE_SLACK_VARIABLES )
        {
            if ( _temporarilyDontUseSlacks )
            {
                log( "Temporarily disabling slacks\n" );
                _temporarilyDontUseSlacks = false;
            }
            else
                computeSlackBounds();
        }

        _reluUpdateFrequency.clear();

        glpkWrapper.setBoundCalculationHook( &boundCalculationHook );
        glpkWrapper.setIterationCountCallback( &iterationCountCallback );
        glpkWrapper.setReportSoiCallback( &reportSoiCallback );

        GlpkWrapper::GlpkAnswer answer;

        try
        {
            answer = glpkWrapper.run( *this );
        }
        catch ( InvariantViolationError &e )
        {
            _currentGlpkWrapper = NULL;
            timeval lpEnd = Time::sampleMicro();
            unsigned timePassed = Time::timePassed( lpStart, lpEnd );
            _totalLpSolverTimeMilli += timePassed;

            if ( timePassed > _maxLpSolverTimeMilli )
                _maxLpSolverTimeMilli = timePassed;

            throw e;
        }

        _currentGlpkWrapper = NULL;
        timeval lpEnd = Time::sampleMicro();
        unsigned timePassed = Time::timePassed( lpStart, lpEnd );
        _totalLpSolverTimeMilli += timePassed;

        if ( timePassed > _maxLpSolverTimeMilli )
            _maxLpSolverTimeMilli = timePassed;

        if ( answer == GlpkWrapper::SOLUTION_FOUND )
        {
            log( "LP solver solved the problem. Updating tableau and assignment\n" );
            ++_numLpSolverFoundSolution;

            timeval extractionStart = Time::sampleMicro();

            // Two options: either restore by basics (and pivot manually), or just take the entire tableau.
            if ( _glpkExtractJustBasics )
            {
                Set<unsigned> newBasics;
                glpkWrapper.extractBasicVariables( *this, newBasics );

                Set<unsigned> shouldBeBasic = Set<unsigned>::difference( newBasics, _basicVariables );
                Set<unsigned> shouldntBeBasic = Set<unsigned>::difference( _basicVariables, newBasics );
                adjustBasicVariables( shouldBeBasic, shouldntBeBasic, false );
            }
            else
            {
                glpkWrapper.extractTableau( this, &_tableau, &_basicVariables, &_eliminatedVars );
            }

            Map<unsigned, double> assignment;
            glpkWrapper.extractAssignment( *this, assignment );
            adjustGlpkAssignment( assignment );

            for ( const auto &pair : assignment )
            {
                double value = pair.second;
                _assignment[pair.first] = value;
            }

            calculateBasicVariableValues();
            computeVariableStatus();

            unsigned brokenReluPairs;
            unsigned brokenNonBasicReluPairs;
            countBrokenReluPairs( brokenReluPairs, brokenNonBasicReluPairs );

            if ( brokenReluPairs > _maxBrokenReluAfterGlpk )
                _maxBrokenReluAfterGlpk = brokenReluPairs;
            _totalBrokenReluAfterGlpk += brokenReluPairs;
            _totalBrokenNonBasicReluAfterGlpk += brokenNonBasicReluPairs;

            timeval extractionEnd = Time::sampleMicro();
            unsigned timePassed = Time::timePassed( extractionStart, extractionEnd );
            _totalLpExtractionTime += timePassed;

            DEBUG( checkInvariants() );

            if ( !allVarsWithinBounds() )
//            if ( !allVarsWithinBounds( true ) )
            {
                // This rarely happens, but when it does - need to restore.
                // I'm guessing this is due to numerical instability when restoring the basics.
                log( "Error! Returned from GLPK but have oob variables\n" );

                ++_numLpSolverIncorrectAssignment;
                restoreTableauFromBackup( _consecutiveGlpkFailureCount < 5 );

                dump();

                if ( _previousGlpkAnswer == GlpkWrapper::SOLVER_FAILED )
                {
                    // Two failures in a row, so a restoration didn't help.
                    // Next time, don't use slacks.
                    _temporarilyDontUseSlacks = true;
                }

                _previousGlpkAnswer = GlpkWrapper::SOLVER_FAILED;
                ++_consecutiveGlpkFailureCount;
                return GlpkWrapper::SOLVER_FAILED;
            }

            _previousGlpkAnswer = GlpkWrapper::SOLUTION_FOUND;
            _consecutiveGlpkFailureCount = 0;
            return GlpkWrapper::SOLUTION_FOUND;
        }
        else if ( answer == GlpkWrapper::NO_SOLUTION_EXISTS )
        {
            log( "LP solver showed no solution exists\n" );
            ++_numLpSolverNoSolution;
            _previousGlpkAnswer = GlpkWrapper::NO_SOLUTION_EXISTS;
            _consecutiveGlpkFailureCount = 0;
            return GlpkWrapper::NO_SOLUTION_EXISTS;
        }

        log( "LP solver failed! Restoring from original matrix...\n" );
        ++_numLpSolverFailed;
        restoreTableauFromBackup( _consecutiveGlpkFailureCount < 5 );

        dump();

        if ( _previousGlpkAnswer == GlpkWrapper::SOLVER_FAILED )
        {
            // Two failures in a row, so a restoration didn't help.
            // Next time, don't use slacks.
            _temporarilyDontUseSlacks = true;
        }

        _previousGlpkAnswer = GlpkWrapper::SOLVER_FAILED;
        ++_consecutiveGlpkFailureCount;
        return GlpkWrapper::SOLVER_FAILED;
    }

    void storeGlpkBoundTightening( int n, int m, int *head,
                                   int leavingBasic,
                                   int enteringNonBasicEncoding,
                                   double *basicRow )
    {
        List<GlpkRowEntry> row;
        row.append( GlpkRowEntry( _currentGlpkWrapper->glpkEncodingToVariable( head[leavingBasic] ), -1.0 ) );

        unsigned numberOfNonBasics = n - m;
        double weightOfEntering = 0.0;
        unsigned enteringNonBasic = _currentGlpkWrapper->glpkEncodingToVariable( head[m + enteringNonBasicEncoding] );

        for ( unsigned i = 1; i <= numberOfNonBasics; ++i )
        {
            unsigned nonBasic = _currentGlpkWrapper->glpkEncodingToVariable( head[i + m] );
            double weight = basicRow[i];

            if ( nonBasic == enteringNonBasic )
                weightOfEntering = weight;

            if ( !FloatUtils::isZero( weight ) )
                row.append( GlpkRowEntry( nonBasic, weight ) );
        }

        if ( !_maximalGlpkBoundTightening )
        {
            storeGlpkBoundTighteningOnRow( row, _currentGlpkWrapper->glpkEncodingToVariable( head[leavingBasic] ) );

            DEBUG({
                      if ( FloatUtils::isZero( weightOfEntering ) )
                      {
                          printf( "Error! weightOfEntering is zero!\n" );
                          exit( 1 );
                      }
                  });

            double scale = -1.0 / weightOfEntering;
            for ( auto &it : row )
            {
                if ( it._variable == enteringNonBasic )
                    it._coefficient = -1.0;
                else
                    it._coefficient *= scale;
            }

            storeGlpkBoundTighteningOnRow( row, enteringNonBasic );
        }
        else
        {
            // If maximal tightening is on, learn all possible bounds from this row
            for ( const auto &newBasic : row )
            {
                List<GlpkRowEntry> copy = row;
                if ( !FloatUtils::areEqual( newBasic._coefficient, -1.0 ) )
                {
                    // Scale the copy
                    double scale = -1.0 / newBasic._coefficient;
                    for ( auto &it : copy )
                    {
                        if ( it._variable == newBasic._variable )
                            it._coefficient = -1.0;
                        else
                            it._coefficient *= scale;
                    }
                }

                storeGlpkBoundTighteningOnRow( copy, newBasic._variable );
            }
        }
    }

    void storeGlpkBoundTighteningOnRow( List<GlpkRowEntry> &row, unsigned basic )
    {
        if ( ( _useSlackVariablesForRelus == USE_ROW_AND_COL_SLACK_VARIABLES ) &&
             _activeSlackRowVars.exists( basic ) )
        {
            // When using row and col slacks, rows are always fixed at 0, so ignore.
            return;
        }

        double max = 0.0;
        double min = 0.0;
        unsigned minBoundLevel = 0;
        unsigned maxBoundLevel = 0;

        ++_storeGlpkBoundTighteningCalls;
        if ( _activeSlackRowVars.exists( basic ) || _activeSlackColVars.exists( basic ) )
            ++_storeGlpkBoundTighteningCallsOnSlacks;

        DEBUG( Set<unsigned> seenVariables; );

        for ( const auto &entry : row )
        {
            if ( entry._variable == basic )
            {
                DEBUG({
                          if ( FloatUtils::areDisequal( entry._coefficient, -1.0 ) )
                          {
                              printf( "Error! storeGlpkBoundTighteningOnRow expected -1.0 coefficient for basic!\n" );
                              exit( 1 );
                          }
                      });

                continue;
            }

            unsigned nonBasic = entry._variable;
            double weight = entry._coefficient;

            // TODO: ignore tiny weights, for numerical stability

            DEBUG({
                      if ( !_activeSlackRowVars.exists( nonBasic ) && !_activeSlackColVars.exists( nonBasic ) &&
                           ( ( _tableau.getColumnSize( nonBasic ) == 0 ) ||
                             _eliminatedVars.exists( nonBasic ) ||
                             isDissolvedBVariable( nonBasic ) ) )
                      {
                          printf( "Error! A non active non-basic variable appeared!\n" );
                          exit( 1 );
                      }

                      if ( seenVariables.exists( nonBasic ) )
                      {
                          printf( "Error! Same variable twice!\n" );
                          exit( 1 );
                      }
                      seenVariables.insert( nonBasic );

                      if ( nonBasic == basic )
                      {
                          printf( "Error: basic == nonbasic!\n" );
                          exit( 1 );
                      }

                      if ( !_activeSlackRowVars.exists( nonBasic ) && !_activeSlackColVars.exists( nonBasic ) &&
                           ( !_lowerBounds[nonBasic].finite() || !_upperBounds[nonBasic].finite() ) )
                      {
                          printf( "Error! Encountered an infinite bound!\n" );
                          exit( 1 );
                      }
                  });

            double currentLowerNB;
            double currentUpperNB;

            unsigned currentLowerNBLevel;
            unsigned currentUpperNBLevel;

            if ( !_activeSlackRowVars.exists( nonBasic ) && !_activeSlackColVars.exists( nonBasic ) )
            {
                if ( _glpkStoredLowerBounds.exists( nonBasic ) )
                {
                    currentLowerNB = _glpkStoredLowerBounds[nonBasic].getBound();
                    currentLowerNBLevel = _glpkStoredLowerBounds[nonBasic].getLevel();
                }
                else
                {
                    currentLowerNB = _lowerBounds[nonBasic].getBound();
                    currentLowerNBLevel = _lowerBounds[nonBasic].getLevel();
                }

                if ( _glpkStoredUpperBounds.exists( nonBasic ) )
                {
                    currentUpperNB = _glpkStoredUpperBounds[nonBasic].getBound();
                    currentUpperNBLevel = _glpkStoredUpperBounds[nonBasic].getLevel();
                }
                else
                {
                    currentUpperNB = _upperBounds[nonBasic].getBound();
                    currentUpperNBLevel = _upperBounds[nonBasic].getLevel();
                }
            }
            else if ( _activeSlackColVars.exists( nonBasic ) )
            {
                currentLowerNB = _slackToLowerBound[nonBasic].getBound();
                currentUpperNB = _slackToUpperBound[nonBasic].getBound();

                currentLowerNBLevel = _slackToLowerBound[nonBasic].getLevel();
                currentUpperNBLevel = _slackToUpperBound[nonBasic].getLevel();
            }
            else
            {
                // Slack row var. If we're using cols, the bounds are 0. Otherwise, read real bounds.
                // This should have no effect on the level of the learned bounds
                if ( _useSlackVariablesForRelus == USE_ROW_AND_COL_SLACK_VARIABLES )
                {
                    currentLowerNB = 0.0;
                    currentUpperNB = 0.0;

                    currentLowerNBLevel = 0.0;
                    currentUpperNBLevel = 0.0;
                }
                else
                {
                    currentLowerNB = _slackToLowerBound[nonBasic].getBound();
                    currentUpperNB = _slackToUpperBound[nonBasic].getBound();

                    currentLowerNBLevel = _slackToLowerBound[nonBasic].getLevel();
                    currentUpperNBLevel = _slackToUpperBound[nonBasic].getLevel();
                }
            }

            if ( FloatUtils::isPositive( weight ) )
            {
                max += ( currentUpperNB * weight );
                min += ( currentLowerNB * weight );

                if ( minBoundLevel < currentLowerNBLevel )
                    minBoundLevel = currentLowerNBLevel;
                if ( maxBoundLevel < currentUpperNBLevel )
                    maxBoundLevel = currentUpperNBLevel;
            }
            else if ( FloatUtils::isNegative( weight ) )
            {
                min += ( currentUpperNB * weight );
                max += ( currentLowerNB * weight );

                if ( maxBoundLevel < currentLowerNBLevel )
                    maxBoundLevel = currentLowerNBLevel;
                if ( minBoundLevel < currentUpperNBLevel )
                    minBoundLevel = currentUpperNBLevel;
            }
        }

        double currentLower;
        double currentUpper;
        double currentLowerLevel;
        double currentUpperLevel;

        if ( _activeSlackColVars.exists( basic ) )
        {
            DEBUG({
                      if ( _useSlackVariablesForRelus != USE_ROW_AND_COL_SLACK_VARIABLES )
                      {
                          printf( "Error! Learned a bound for a col slack variable!\n" );
                          exit( 1 );
                      }
                  });

            currentLower = _slackToLowerBound[basic].getBound();
            currentLowerLevel = _slackToLowerBound[basic].getLevel();
            currentUpper = _slackToUpperBound[basic].getBound();
            currentUpperLevel = _slackToUpperBound[basic].getLevel();
        }
        else if ( _activeSlackRowVars.exists( basic ) )
        {
            DEBUG({
                      if ( _useSlackVariablesForRelus != USE_ROW_SLACK_VARIABLES )
                      {
                          printf( "Error! Learned a bound for a row slack variable!\n" );
                          exit( 1 );
                      }
                  });

            currentLower = _slackToLowerBound[basic].getBound();
            currentLowerLevel = _slackToLowerBound[basic].getLevel();
            currentUpper = _slackToUpperBound[basic].getBound();
            currentUpperLevel = _slackToUpperBound[basic].getLevel();
        }
        else
        {
            if ( _glpkStoredLowerBounds.exists( basic ) )
            {
                currentLower = _glpkStoredLowerBounds[basic].getBound();
                currentLowerLevel = _glpkStoredLowerBounds[basic].getLevel();
            }
            else
            {
                currentLower = _lowerBounds[basic].getBound();
                currentLowerLevel = _lowerBounds[basic].getLevel();
            }

            if ( _glpkStoredUpperBounds.exists( basic ) )
            {
                currentUpper = _glpkStoredUpperBounds[basic].getBound();
                currentUpperLevel = _glpkStoredUpperBounds[basic].getLevel();
            }
            else
            {
                currentUpper = _upperBounds[basic].getBound();
                currentUpperLevel = _upperBounds[basic].getLevel();
            }
        }

        bool updateOccurred = false;
        if ( FloatUtils::lt( max, currentUpper ) )
        {
            if ( _activeSlackColVars.exists( basic ) || _activeSlackRowVars.exists( basic ) )
            {
                _slackToUpperBound[basic].setBound( max );
                _slackToUpperBound[basic].setLevel( maxBoundLevel );
            }
            else
            {
                _glpkStoredUpperBounds[basic].setBound( max );
                _glpkStoredUpperBounds[basic].setLevel( maxBoundLevel );
            }

            ++_numBoundsDerivedThroughGlpk;
            updateOccurred = true;
            currentUpper = max;
            currentUpperLevel = maxBoundLevel;
        }

        if ( FloatUtils::gt( min, currentLower ) )
        {
            if ( _activeSlackColVars.exists( basic ) || _activeSlackRowVars.exists( basic ) )
            {
                _slackToLowerBound[basic].setBound( min );
                _slackToLowerBound[basic].setLevel( minBoundLevel );
            }
            else
            {
                _glpkStoredLowerBounds[basic].setBound( min );
                _glpkStoredLowerBounds[basic].setLevel( minBoundLevel );
            }

            ++_numBoundsDerivedThroughGlpk;
            updateOccurred = true;
            currentLower = min;
            currentLowerLevel = minBoundLevel;
        }

        if ( updateOccurred && FloatUtils::gt( currentLower, currentUpper ) )
        {
            throw InvariantViolationError( std::max( currentLowerLevel, currentUpperLevel ) );
        }
    }

    bool learnedGlpkBounds() const
    {
        return _glpkStoredLowerBounds.size() + _glpkStoredUpperBounds.size() > 0;
    }

    void performGlpkBoundTightening()
    {
        log( "Starting GLPK bound tightening\n" );

        // It is wrong to assume that all bounds are improvements over existing bounds. As
        // we begin updating, things may change because of relu stuff - so check again before
        // each bound.
        for ( const auto &lowerBound : _glpkStoredLowerBounds )
        {
            if ( FloatUtils::gt( lowerBound.second.getBound(), _lowerBounds[lowerBound.first].getBound() ) )
                updateLowerBound( lowerBound.first, lowerBound.second.getBound(), lowerBound.second.getLevel() );
        }

        for ( const auto &upperBound : _glpkStoredUpperBounds )
        {
            if ( FloatUtils::lt( upperBound.second.getBound(), _upperBounds[upperBound.first].getBound() ) )
                updateUpperBound( upperBound.first, upperBound.second.getBound(), upperBound.second.getLevel() );
        }

        _glpkStoredLowerBounds.clear();
        _glpkStoredUpperBounds.clear();

        log( "Finished with GLPK bound tightening\n" );
    }

    void glpkIterationCountCallback( int count )
    {
        log( Stringf( "GLPK: number of iterations = %i\n", count ) );
        _totalLpPivots += count;
    }

    void glpkReportSoi( double soi )
    {
        log( Stringf( "GLPK report soi: %.10lf\n", soi ) );
        _glpkSoi = soi;
    }

    void setUpperBound( unsigned variable, double bound )
    {
        _upperBounds[variable].setBound( bound );
        _upperBounds[variable].setLevel( 0 );
    }

    void setLowerBound( unsigned variable, double bound )
    {
        _lowerBounds[variable].setBound( bound );
        _lowerBounds[variable].setLevel( 0 );
    }

    bool boundInvariantHolds( unsigned variable, unsigned &violatingStackLevel )
    {
        //判断已设定的上下界值是否是合理的，即lowerBounds要小于等于upperBounds，
        //如果其中一边为无穷大，由于无穷大是无法比较的，则默认都为true

        if ( !_upperBounds[variable].finite() || !_lowerBounds[variable].finite() )
            return true;

        // 如果不是，并返回false，由调用函数抛出异常
        if ( !FloatUtils::lte( _lowerBounds[variable].getBound(), _upperBounds[variable].getBound() ) )
        {
            violatingStackLevel = std::max( _lowerBounds[variable].getLevel(), _upperBounds[variable].getLevel() );
            return false;
        }

        // 其他所有情况都返回true
        return true;
    }

    bool activeReluVariable( unsigned variable ) const
    {
        if ( !_reluPairs.isRelu( variable ) )
            return false;

        unsigned f = _reluPairs.isF( variable ) ? variable : _reluPairs.toPartner( variable );
        return !_dissolvedReluVariables.exists( f );
    }

    void updateUpperBound( unsigned variable, double bound, unsigned level )
    {
//        printf( "\nenter updateUpperBound: variable:%u, bound:%.10f\n", variable, bound );

        unsigned partner = 0, b = 0, f = 0;

        // 如果是relu变量，则取出它以及相关的b或f
        if ( _reluPairs.isRelu( variable ) )
        {
            // The variable is relu.
            partner = _reluPairs.toPartner( variable );
            f = _reluPairs.isF( variable ) ? variable : partner;
            b = _reluPairs.isB( variable ) ? variable : partner;
        }

        // 如果不是relu变量，或者relu变量的f已经被消除
        if ( !_reluPairs.isRelu( variable ) || _dissolvedReluVariables.exists( f ) )
        {
            // For non-relus, we can just update the bound.
            // 直接更新上界为传入的bound值
            _upperBounds[variable].setBound( bound );
            _upperBounds[variable].setLevel( level );

            unsigned violatingStackLevel;

            // 判断上下界是否是符合常理的
            if ( !boundInvariantHolds( variable, violatingStackLevel ) )
                throw InvariantViolationError( violatingStackLevel );

            computeVariableStatus( variable );

            // If the variable is basic, it's okay if it's out of bounds.
            // If non-basic and out of bounds, need to update.
            // 如果更新了界限范围的是一个non-basic，且越界了，那么就要将其值修正为这个界限值（由于update可以接受的delta可为正负数，则这里可以直接用bound - _assignment[variable]
            if ( !_basicVariables.exists( variable ) && outOfBounds( variable ) )
                update( variable, bound - _assignment[variable] );

            return;
        }

        // We are in the relu case. Should we use almost-broken elimination?
        if ( FloatUtils::isPositive( bound ) &&
             FloatUtils::lte( bound, ALMOST_BROKEN_RELU_MARGIN ) )
        {
            ++_almostBrokenReluPairCount;

            if ( _eliminateAlmostBrokenRelus )
            {
                ++_almostBrokenReluPairFixedCount;
                bound = 0.0;
            }
        }

        // 如果上界是正数，更新上界，如果越界，也是更新为正数
        // If the bound is positive, update bounds on both F and B.
        if ( FloatUtils::isPositive( bound ) )
        {
            _upperBounds[variable].setBound( bound );
            _upperBounds[variable].setLevel( level );
            _upperBounds[partner].setBound( bound );
            _upperBounds[partner].setLevel( level );

            unsigned violatingStackLevel;
            if ( !boundInvariantHolds( variable, violatingStackLevel ) ||
                 !boundInvariantHolds( partner, violatingStackLevel ) )
                throw InvariantViolationError( violatingStackLevel );

            computeVariableStatus( variable );
            computeVariableStatus( partner );

            // Violations are okay for basic, but need to update if non-basic
            if ( !_basicVariables.exists( variable ) && outOfBounds( variable ) )
                update( variable, bound - _assignment[variable], true );
            if ( !_basicVariables.exists( partner ) && outOfBounds( partner ) )
                update( partner, bound - _assignment[partner], true );

            return;
        }
        else
        {
//            printf("\nDo a split~~~~\n");

            // 如果上界是0或负数，就可以假设已经求解了这个f，并设置
            // 此时如果越界，也是更新为上界（负数），那么就是要乘以ratio
//            printf( "enter updateUpperBound---- upperBound is zero or negative\n" );

//            _upperBounds[variable].setBound( bound );
//            _upperBounds[variable].setLevel( level );
//            _upperBounds[partner].setBound( bound );
//            _upperBounds[partner].setLevel( level );

            // 仅仅是更新界限
            if (_reluPairs.isF(variable)) {     // 如果要更新的界限是f的，那么b就是 （ 1 / leakyRatio）* bound
                _upperBounds[variable].setBound( bound );
                _upperBounds[variable].setLevel( level );
                _upperBounds[partner].setBound( (1 / getLeakyValue()) * bound );
                _upperBounds[partner].setLevel( level );
            } else{ // 如果要更新的界限是b的，那么f的就是 leakyRatio * bound
                _upperBounds[variable].setBound( bound );
                _upperBounds[variable].setLevel( level );
                _upperBounds[partner].setBound( getLeakyValue() * bound );
                _upperBounds[partner].setLevel( level );
            }

            /*****  上界为负数，则取值一定是负数，既然要mark为resolved，那么就要连同assignment一起改变，还要添加新的约束 *****/

            // 将要更新的值设置为非基变量，因为如果是基变量，就不会更新联动的值，第二个参数表示禁止将其设置为pivot的候选对象
            if ( _basicVariables.exists( b ) )
            {
//                printf("~~~~~Before add split restrain, the b is basic , so need to change to non-basic\n");
                makeNonBasic( b, f );
            }
            if ( _basicVariables.exists( f ) )
            {
//                printf("~~~~~Before add split restrain, the f is basic , so need to change to non-basic\n");
                makeNonBasic( f, b );
            }
            dump();
            // b和f在整个运算过程中应该都存在，除非是dissolve，此时如果出现这种情况，是进入不到这里来的，但如果存在，就报错
            if (_tableau.activeColumn(b) && _tableau.activeColumn(f)) {
                _tableau.eraseRow(f);
                _tableau.replaceNonBasicWithAnotherNonBasic(f, b, getLeakyValue());
//                deleteIfJustItself();  //上面一行有过tableau的合并，可能会导致某些var的行只有它自己，跟其他变量不再有关系，所以要检测，删除这样的行
                _tableau.addEntry(f, b, getLeakyValue());
                _tableau.addEntry(f, f, -1);
                markBasic(f);
            } else{
                throw Error( Error::DISSOLVED_SPLIT_ERROR );
            }
//            printf("~~~~~After add new restrain");
            dump();
            // 将b设置为与f匹配的值，如果越界，后面会处理
            update( b, ( 1 /getLeakyValue() ) * _assignment[f] - _assignment[b], true );

            /***** 修复可能出现的错误值 ****/
            // 更新完之后，与b相关联的行的赋值，可能会出现错误，此时要以b为标准，进行值的修复
            for (unsigned row = 0; row < _tableau.getNumVars(); row++) {
                // 如果系数不为0，则b在这一行有出现过，重新计算这一行的值
//                if ( (row != f) && (!FloatUtils::areEqual(_tableau.getCell(row, b), 0.0)) ){
                // 注意：这里重新赋值是包括f那一行，因为在update b与f一致的过程中，可能又因为牵连会导致f改变，为了确保计算一致，这里也要更改
                if (_tableau.getCell(row, b)!= 0.0 ){
                    double sum = 0.0;
                    for (unsigned col = 0; col < _tableau.getNumVars(); col++) {
                        // 如果有值，那就取出计算,
                        // 一般col和row相同的时候用-1表示自身，所以要跳过col==row的情况
//                        if ((col != row) && (!FloatUtils::areEqual(_tableau.getCell(row, col), 0.0))) {
                        if ((row != col) && (_tableau.getCell(row, col)!= 0.0)) {
                            sum += _tableau.getCell(row, col) * _assignment[col];
                        }
                    }
                    _assignment[row] = sum;
                }
            }

//            printf("~~~~~After update b to f\n");
            dump();

            markReluVariableDissolved( f, TYPE_SPLIT );

//            printf( "~~~~~_assignment[b]: %.10f, _assignment[f]: %.10f\n", _assignment[b], _assignment[f] );

            /*** 其他设置 ****/
            unsigned violatingStackLevel;
            if ( !boundInvariantHolds( variable, violatingStackLevel ) ||
                 !boundInvariantHolds( partner, violatingStackLevel ) )
                throw InvariantViolationError( violatingStackLevel );

            computeVariableStatus( variable );
            computeVariableStatus( partner );

            // Violations are okay for basic, but need to update if non-basic
            if ( !_basicVariables.exists( variable ) && outOfBounds( variable ) ){
//                printf("~~~~~After split, some variable is out of bound\n");
                update( variable, bound - _assignment[variable], true );
            }
            if ( !_basicVariables.exists( partner ) && outOfBounds( partner ) ){
//                printf("~~~~~After split, some variable is out of bound\n");
                update( partner, bound - _assignment[partner], true );
            }

            // Violations are okay for basic, but need to update if non-basic
//            if ( !_basicVariables.exists( b ) && outOfBounds( b ) )
//                update( b, bound - _assignment[b], true );
//            if ( !_basicVariables.exists( f ) && outOfBounds( f ) )
//                update( f, leakyRatio * bound - _assignment[f], true );

            return;
        }
    }

    void updateUpperBound_temp( unsigned variable, double bound, unsigned level )
    {
        unsigned partner = 0, b = 0, f = 0;

        // 如果是relu变量，则取出它以及相关的b或f
        if ( _reluPairs.isRelu( variable ) )
        {
            // The variable is relu.
            partner = _reluPairs.toPartner( variable );
            f = _reluPairs.isF( variable ) ? variable : partner;
            b = _reluPairs.isB( variable ) ? variable : partner;
        }

        // 如果不是relu变量，或者relu变量的f已经被消除
        if ( !_reluPairs.isRelu( variable ) || _dissolvedReluVariables.exists( f ) )
        {
            // For non-relus, we can just update the bound.
            // 直接更新上界为传入的bound值
            _upperBounds[variable].setBound( bound );
            _upperBounds[variable].setLevel( level );

            unsigned violatingStackLevel;

            // 判断上下界是否是符合常理的
            if ( !boundInvariantHolds( variable, violatingStackLevel ) )
                throw InvariantViolationError( violatingStackLevel );

            computeVariableStatus( variable );

            // If the variable is basic, it's okay if it's out of bounds.
            // If non-basic and out of bounds, need to update.
            // 如果更新了界限范围的是一个non-basic，且越界了，那么就要将其值修正为这个界限值（由于update可以接受的delta可为正负数，则这里可以直接用bound - _assignment[variable]
            if ( !_basicVariables.exists( variable ) && outOfBounds( variable ) )
                update( variable, bound - _assignment[variable] );

            return;
        }

        // We are in the relu case. Should we use almost-broken elimination?
        if ( FloatUtils::isPositive( bound ) &&
             FloatUtils::lte( bound, ALMOST_BROKEN_RELU_MARGIN ) )
        {
            ++_almostBrokenReluPairCount;

            if ( _eliminateAlmostBrokenRelus )
            {
                ++_almostBrokenReluPairFixedCount;
                bound = 0.0;
            }
        }

        // If the bound is positive, update bounds on both F and B.
        if ( FloatUtils::isPositive( bound ) )
        {
            _upperBounds[variable].setBound( bound );
            _upperBounds[variable].setLevel( level );
            _upperBounds[partner].setBound( bound );
            _upperBounds[partner].setLevel( level );

            unsigned violatingStackLevel;
            if ( !boundInvariantHolds( variable, violatingStackLevel ) ||
                 !boundInvariantHolds( partner, violatingStackLevel ) )
                throw InvariantViolationError( violatingStackLevel );

            computeVariableStatus( variable );
            computeVariableStatus( partner );

            // Violations are okay for basic, but need to update if non-basic
            if ( !_basicVariables.exists( variable ) && outOfBounds( variable ) )
                update( variable, bound - _assignment[variable], true );
            if ( !_basicVariables.exists( partner ) && outOfBounds( partner ) )
                update( partner, bound - _assignment[partner], true );

            return;
        }
        else
        {
            // Non-positive bound.
            if ( FloatUtils::isNegative( bound ) && _reluPairs.isF( variable ) )
            {
                _upperBounds[variable].setBound( bound );
                _upperBounds[variable].setLevel( level );

                // F variables have lower bound of at least 0. Do this "by-the-book" for proper stack handling.
                unsigned violatingStackLevel;
                if ( !boundInvariantHolds( variable, violatingStackLevel ) )
                    throw InvariantViolationError( violatingStackLevel );
                else
                {
                    printf( "Error! Expected violation on F!\n" );
                    exit( 1 );
                }
            }

            // F will have zero as upper bound.
            markReluVariableDissolved( f, TYPE_SPLIT );

            _upperBounds[f].setBound( 0.0 );
            _upperBounds[f].setLevel( level );
            _upperBounds[b].setBound( bound );
            _upperBounds[b].setLevel( level );

            unsigned violatingStackLevel;
            if ( !boundInvariantHolds( b, violatingStackLevel ) || !boundInvariantHolds( f, violatingStackLevel ) )
                throw InvariantViolationError( violatingStackLevel );

            computeVariableStatus( b );
            computeVariableStatus( f );

            // Violations are okay for basic, but need to update if non-basic
            if ( !_basicVariables.exists( b ) && outOfBounds( b ) )
                update( b, bound - _assignment[b], true );
            if ( !_basicVariables.exists( f ) && outOfBounds( f ) )
                update( f, -_assignment[f], true );

            return;
        }
    }

    bool updateLowerBound( unsigned variable, double bound, unsigned level )
    {
//        printf( "\nenter updateLowerBound: variable:%u, bound:%.10f\n", variable, bound );

        unsigned partner = 0, b = 0, f = 0;

        // 如果是relu变量，找到b\f
        if ( _reluPairs.isRelu( variable ) )
        {
            // The variable is relu.
            partner = _reluPairs.toPartner( variable );
            f = _reluPairs.isF( variable ) ? variable : partner;
            b = _reluPairs.isB( variable ) ? variable : partner;
        }

        // 如果不是relu变量，更新相应的下界
        if ( !_reluPairs.isRelu( variable ) || _dissolvedReluVariables.exists( f ) )
        {
            // For non-relus, we can just update the bound.
            _lowerBounds[variable].setBound( bound );
            _lowerBounds[variable].setLevel( level );

            unsigned violatingStackLevel;
            if ( !boundInvariantHolds( variable, violatingStackLevel ) )
                throw InvariantViolationError( violatingStackLevel );

            computeVariableStatus( variable );

            // If the variable is basic, it's okay if it's out of bounds.
            // If non-basic, need to update.
            if ( !_basicVariables.exists( variable ) && outOfBounds( variable ) )
                update( variable, bound - _assignment[variable] );

            // The tableau has not changed
            return false;
        }

        // 没开启，忽略
        // Should we use almost-broken elimination?
        if ( FloatUtils::isNegative( bound ) &&
             FloatUtils::gte( bound, -ALMOST_BROKEN_RELU_MARGIN ) )
        {
            ++_almostBrokenReluPairCount;

            if ( _eliminateAlmostBrokenRelus )
            {
                ++_almostBrokenReluPairFixedCount;
                bound = 0.0;
            }
        }

        // 如果是relu变量，且下界大于等于0，则更新b与f的下界，并且可以尝试设置求解了f为merge
        // 此时即使越界，要更新也是更新为正数，所以update时只是到bound
        // If the bound is non-negative, update bounds on both F and B.
        if ( !FloatUtils::isNegative( bound ) )
        {
            log( "Update lower bound: non-negative lower bound\n" );

//            printf("\nDo a merge~~~~\n");

            _lowerBounds[variable].setBound( bound );
            _lowerBounds[variable].setLevel( level );
            _lowerBounds[partner].setBound( bound );
            _lowerBounds[partner].setLevel( level );

            unsigned violatingStackLevel;
            if ( !boundInvariantHolds( variable, violatingStackLevel ) ||
                 !boundInvariantHolds( partner, violatingStackLevel ) )
                throw InvariantViolationError( violatingStackLevel );

            computeVariableStatus( variable );
            computeVariableStatus( partner );

            // Violations are okay for basic, but need to update if non-basic
            if ( !_basicVariables.exists( variable ) && outOfBounds( variable ) )
                update( variable, bound - _assignment[variable], true );
            if ( !_basicVariables.exists( partner ) && outOfBounds( partner ) )
                update( partner, bound - _assignment[partner], true );

            // 这里会将b设置为f的值，并将f设置为resolved,如果上面更新的界限导致basic越界，也会在里面进行pivot，然后更新回界限，再合并
            return unifyReluPair( f );
        }
        else
        {
            // 下界为负数时，也还是同时更新两者的下界,此时即使越界，也是更新回下界（负数），所以要乘ratio

            if (_reluPairs.isF(variable)) {     // 如果要更新的界限是f的，那么b就是 （ 1 / leakyRatio）* bound
                _upperBounds[variable].setBound( bound );
                _upperBounds[variable].setLevel( level );
                _upperBounds[partner].setBound( (1 / getLeakyValue()) * bound );
                _upperBounds[partner].setLevel( level );
            } else{ // 如果要更新的界限是b的，那么f的就是 leakyRatio * bound
                _upperBounds[variable].setBound( bound );
                _upperBounds[variable].setLevel( level );
                _upperBounds[partner].setBound( getLeakyValue() * bound );
                _upperBounds[partner].setLevel( level );
            }

            unsigned violatingStackLevel;
            if ( !boundInvariantHolds( variable, violatingStackLevel ) ||
                 !boundInvariantHolds( partner, violatingStackLevel ) )
                throw InvariantViolationError( violatingStackLevel );

            computeVariableStatus( variable );
            computeVariableStatus( partner );


            // Violations are okay for basic, but need to update if non-basic
            if ( !_basicVariables.exists( variable ) && outOfBounds( variable ) )
                update( variable, bound - _assignment[variable], true );
            if ( !_basicVariables.exists( partner ) && outOfBounds( partner ) )
                update( partner, bound - _assignment[partner], true );

            // The tableau has not changed
            return false;
        }
    }
    bool updateLowerBound_temp( unsigned variable, double bound, unsigned level )
    {
        unsigned partner = 0, b = 0, f = 0;

        // 如果是relu变量，找到b\f
        if ( _reluPairs.isRelu( variable ) )
        {
            // The variable is relu.
            partner = _reluPairs.toPartner( variable );
            f = _reluPairs.isF( variable ) ? variable : partner;
            b = _reluPairs.isB( variable ) ? variable : partner;
        }

        // 如果不是relu变量，更新相应的下界
        if ( !_reluPairs.isRelu( variable ) || _dissolvedReluVariables.exists( f ) )
        {
            // For non-relus, we can just update the bound.
            _lowerBounds[variable].setBound( bound );
            _lowerBounds[variable].setLevel( level );

            unsigned violatingStackLevel;
            if ( !boundInvariantHolds( variable, violatingStackLevel ) )
                throw InvariantViolationError( violatingStackLevel );

            computeVariableStatus( variable );

            // If the variable is basic, it's okay if it's out of bounds.
            // If non-basic, need to update.
            if ( !_basicVariables.exists( variable ) && outOfBounds( variable ) )
                update( variable, bound - _assignment[variable] );

            // The tableau has not changed
            return false;
        }

        // 没开启，忽略
        // Should we use almost-broken elimination?
        if ( FloatUtils::isNegative( bound ) &&
             FloatUtils::gte( bound, -ALMOST_BROKEN_RELU_MARGIN ) )
        {
            ++_almostBrokenReluPairCount;

            if ( _eliminateAlmostBrokenRelus )
            {
                ++_almostBrokenReluPairFixedCount;
                bound = 0.0;
            }
        }

        // 如果是relu变量，且是非负数下界，则更新b与f的下界
        // If the bound is non-negative, update bounds on both F and B.
        if ( !FloatUtils::isNegative( bound ) )
        {
//            log( "Update lower bound: non-negative lower bound\n" );

            _lowerBounds[variable].setBound( bound );
            _lowerBounds[variable].setLevel( level );
            _lowerBounds[partner].setBound( bound );
            _lowerBounds[partner].setLevel( level );

            unsigned violatingStackLevel;
            if ( !boundInvariantHolds( variable, violatingStackLevel ) ||
                 !boundInvariantHolds( partner, violatingStackLevel ) )
                throw InvariantViolationError( violatingStackLevel );

            computeVariableStatus( variable );
            computeVariableStatus( partner );

            // Violations are okay for basic, but need to update if non-basic
            if ( !_basicVariables.exists( variable ) && outOfBounds( variable ) )
                update( variable, bound - _assignment[variable], true );
            if ( !_basicVariables.exists( partner ) && outOfBounds( partner ) )
                update( partner, bound - _assignment[partner], true );

            return unifyReluPair( f );
        }
        else
        {
            // Negative bound. This can only be called for the B, doesn't affect the F.

            // 如果是relu变量，且下界是负数，则只更新b的下界，不影响f的下界
            _lowerBounds[variable].setBound( bound );
            _lowerBounds[variable].setLevel( level );

            unsigned violatingStackLevel;
            if ( !boundInvariantHolds( variable, violatingStackLevel ) )
                throw InvariantViolationError( violatingStackLevel );

            computeVariableStatus( variable );

            // If the variable is basic, it's okay if it's out of bounds.
            // If non-basic, need to update.
            if ( !_basicVariables.exists( variable ) && outOfBounds( variable ) )
                update( variable, bound - _assignment[variable], true );

            // The tableau has not changed
            return false;
        }
    }

    // Return true iff the tableau changes
    bool unifyReluPair( unsigned f )
    {
        unsigned b = _reluPairs.toPartner( f );
        log( Stringf( "UnifyReluPair called with f = %s, b = %s\n", toName( f ).ascii(),
                      toName( b ).ascii()) );

        // If these two have been unified before, b's column will be empty. Tableau doesn't change.
        if ( _tableau.getColumnSize( b ) == 0 )
        {
            log( Stringf( "UnifyReluPair: b's column is empty, ignroing. Previous dissolved? %s\n",
                          _dissolvedReluVariables.exists( f ) ? "YES" : "NO" ) );
            return false;
        }

        log( Stringf( "Unifying relu pair: %s, %s\n", toName( b ).ascii(), toName( f ).ascii() ) );

        // First step: make sure f and b are not basic.
        // Note: this may temporarily break the axiom that non-basic variables must be within bounds
        //       (if f or b are currently OOB). However, this will be fixed afterwards.

        if ( _basicVariables.exists( b ) )
            makeNonBasic( b, f );

        if ( _basicVariables.exists( f ) )
            makeNonBasic( f, b );

        log( "Both variables are now non-basic\n" );
        dump();

//        printf("update out of Bounds variables, and set b = f\n");
        // Next: set f to be in bounds.
        if ( tooLow( f ) )
            update( f, _lowerBounds[f].getBound() - _assignment[f], true );
        else if ( tooHigh( f ) )
            update( f, _upperBounds[f].getBound() - _assignment[f], true );

        // Get b to equal f (bounds are equal for both, so this is okay)
        update( b, _assignment[f] - _assignment[b], true );

        // b and f are now equal and non basic. Replace b with f
        _tableau.addColumnEraseSource( b, f );

        // 兼具向_dissolvedReluVariables中添加的功能，
        markReluVariableDissolved( f, TYPE_MERGE );
//        deleteIfJustItself();  //上面一行有过tableau的合并，可能会导致某些var的行只有它自己，跟其他变量不再有关系，所以要检测，删除这样的行

        log( "Tableau after unification:\n" );
        dump();

        return true;
    }

    bool checkIfJustItself(){
        // 首先，只有basic会在对角线上有-1，所以其实只需要遍历basic所在的行
        Set<unsigned> basicSet = getBasicVariables();
        Set<unsigned>::const_iterator itr;
        for(itr = basicSet.begin(); itr != basicSet.end(); itr++){
            bool eraseFlag = true;
            // 遍历basic的一行
            for (unsigned col = 0; col < _tableau.getNumVars(); col++) {
                if ((*itr != col) && (  !FloatUtils::isZero(_tableau.getCell(*itr, col)) ) ){   //如果除了对角线上的地方还有其他地方有元素，
                    // 如果这行中，有除了对角线元素以外的其他系数，就说明是正常的，跳出检测下一行
                    eraseFlag = false;
                    break;
                } else if ((*itr == col) && (_tableau.getCell(*itr, col) != -1) ){  // 如果对角线上不等于-1，出错
                    throw Error( Error::VARIABLE_NOT_BASIC );
                }
            }
            if (eraseFlag) {    // 刚刚检查的一行只有对角线上一个元素为-1，所以是有问题的
//                printf("~~~~~checkIfJustItself(): the row %u is abnormal %s:%u\n", *itr, toName(*itr).ascii(), *itr);
                return true;
            }
        }
        // 检测完全部，发现都是正常的，就返回false
        return false;
    }


    // add by lzs
    void deleteIfJustItself(){
//        printf("\n~~~~~call deleteIfJustItself()\n");
        // 首先，只有basic会在对角线上有-1，所以其实只需要遍历basic所在的行
        Set<unsigned> basicSet = getBasicVariables();
        Set<unsigned>::const_iterator itr;
        for(itr = basicSet.begin(); itr != basicSet.end(); itr++){
            bool eraseFlag = true;
            for (unsigned col = 0; col < _tableau.getNumVars(); col++) {
                if ((*itr != col) && (  !FloatUtils::isZero(_tableau.getCell(*itr, col)) ) ){   //如果除了对角线上的地方还有其他地方有元素，
                    eraseFlag = false;
                    break;
                } else if ((*itr == col) && (_tableau.getCell(*itr, col) != -1) ){  // 如果对角线上不等于-1，出错
                    throw Error( Error::VARIABLE_NOT_BASIC );
                }
            }
            if (eraseFlag) {    // 只有对角线上一个元素为-1，所以要删除
//                printf("~~~~~ delete %s:%u\n", toName(*itr).ascii(), *itr);
                _tableau.eraseRow(*itr);
                _basicVariables.erase(*itr);
            }
        }
    }

    // Return true iff the tableau changes
    bool unifyReluPair_temp( unsigned f )
    {
        unsigned b = _reluPairs.toPartner( f );
        log( Stringf( "UnifyReluPair called with f = %s, b = %s\n", toName( f ).ascii(),
                      toName( b ).ascii()) );

        // If these two have been unified before, b's column will be empty. Tableau doesn't change.
        if ( _tableau.getColumnSize( b ) == 0 )
        {
            log( Stringf( "UnifyReluPair: b's column is empty, ignroing. Previous dissolved? %s\n",
                          _dissolvedReluVariables.exists( f ) ? "YES" : "NO" ) );
            return false;
        }

        log( Stringf( "Unifying relu pair: %s, %s\n", toName( b ).ascii(), toName( f ).ascii() ) );

        // First step: make sure f and b are not basic.
        // Note: this may temporarily break the axiom that non-basic variables must be within bounds
        //       (if f or b are currently OOB). However, this will be fixed afterwards.

        if ( _basicVariables.exists( b ) )
            makeNonBasic( b, f );

        if ( _basicVariables.exists( f ) )
            makeNonBasic( f, b );

        log( "Both variables are now non-basic\n" );
        dump();

        // Next: set f to be in bounds.
        if ( tooLow( f ) )
            update( f, _lowerBounds[f].getBound() - _assignment[f], true );
        else if ( tooHigh( f ) )
            update( f, _upperBounds[f].getBound() - _assignment[f], true );

        // Get b to equal f (bounds are equal for both, so this is okay)
        update( b, _assignment[f] - _assignment[b], true );

        // b and f are now equal and non basic. Replace b with f
        _tableau.addColumnEraseSource( b, f );

        // 兼具向_dissolvedReluVariables中添加的功能，
        markReluVariableDissolved( f, TYPE_MERGE );

        log( "Tableau after unification:\n" );
        dump();

        return true;
    }

    void makeNonBasic( unsigned basic, unsigned forbiddenPartner )
    {
//        printf("enter make non-basic for basic: %u\n", basic);
        if ( !_basicVariables.exists( basic ) )
            throw Error( Error::VARIABLE_NOT_BASIC );

        const Tableau::Entry *rowEntry = _tableau.getRow( basic );
        const Tableau::Entry *current;

        bool found = false;
        unsigned leastEvilNonBasic = 0;
        double leastEvilWeight = 0.0;

        while ( rowEntry )
        {
            current = rowEntry;
            rowEntry = rowEntry->nextInRow();

            unsigned column = current->getColumn();
            if ( ( column == basic ) || ( column == forbiddenPartner ) )
                continue;

            double weight = FloatUtils::abs( current->getValue() );
            if ( FloatUtils::gte( weight, NUMBERICAL_INSTABILITY_CONSTANT ) )
            {
                pivot( column, basic );
                return;
            }

            // The weight is smaller than the instability constant
            found = true;
            if ( FloatUtils::gt( weight, leastEvilWeight ) )
            {
                leastEvilWeight = weight;
                leastEvilNonBasic = column;
            }
        }

        if ( !found )
            throw Error( Error::CANT_MAKE_NON_BASIC );

        pivot( leastEvilNonBasic, basic );
    }

    void setReluPair( unsigned backward, unsigned forward )
    {
        _reluPairs.addPair( backward, forward );

        if ( _useSlackVariablesForRelus != DONT_USE_SLACK_VARIABLES )
        {
            unsigned nextIndex = _fToSlackRowVar.size() + _fToSlackColVar.size();
            _fToSlackRowVar[forward] = _numVariables + nextIndex;
            _slackRowVariableToF[_numVariables + nextIndex] = forward;
            _slackRowVariableToB[_numVariables + nextIndex] = backward;

            if ( _useSlackVariablesForRelus == USE_ROW_AND_COL_SLACK_VARIABLES )
            {
                ++nextIndex;
                _fToSlackColVar[forward] = _numVariables + nextIndex;
            }
        }
    }

    void initializeCell( unsigned row, unsigned column, double value )
    {
        _tableau.addEntry( row, column, value );
    }

    void markBasic( unsigned variable )
    {
        _basicVariables.insert( variable );
    }

    void setName( unsigned variable, String name )
    {
        log( Stringf( "Setting name: %s --> %u\n", name.ascii(), variable ) );
        _variableNames[variable] = name;
    }

    String toName( unsigned variable ) const
    {
        if ( _variableNames.exists( variable ) )
            return _variableNames.at( variable );
        return Stringf( "%u", variable );
    }

    void initialUpdate()
    {
        // 计算当前变量值是否越界的状态,write in _varToStatus, so we can judge their bounds more convenient
        computeVariableStatus();

        for ( unsigned i = 0; i < _numVariables; ++i )
        {
            unsigned violatingStackLevel;

            // 判断所有值的上下界是否是符合常理的(That is : lower bound is less than upper bound)，如果不，要抛出异常
            if ( !boundInvariantHolds( i, violatingStackLevel ) )
            {
                printf( "Bound invariant violation on variable: %s\n", toName( i ).ascii() );
                printf( "Lower bound = %.5lf, upper bound = %.5lf\n",
                        _lowerBounds[i].getBound(), _upperBounds[i].getBound() );
                throw InvariantViolationError( violatingStackLevel );
            }

            // 如果non-basic变量超出上或下界（OutOfBounds会同时检查上下界，任何一方越界都会返回true），那么要进行更新，
            // 超越下界就更新为下界值，否则更新为上界值。

            // 此处是 paper第一步，v31越界，在此处被更新为下界值

            if ( !_basicVariables.exists( i ) && outOfBounds( i ) )
            {
                // 如果是越下界（严格小于），那就将赋值加上差值，使赋值为下界值
                if ( tooLow( i ) )
                    update( i, _lowerBounds[i].getBound() - _assignment[i] );
                    // 如果是越上界，则赋值为上界值
                else
                    update( i, _upperBounds[i].getBound() - _assignment[i] );
            }
        }

        log( "Checking invariants after initial update\n" );
        // check Tableau的值是否有异常、non-basic是否在上下界范围内、relu是否正常
        checkInvariants();  //if has error, it will exit(1) immediately
    }

    // 先设置好fix时不同variable需要的参数，然后进行fix，优先fix b与f一致，其次再考虑设置f与b一致
    bool fixBrokenRelu( unsigned toFix )
    {
        bool isF = _reluPairs.isF( toFix );
        unsigned partner = _reluPairs.toPartner( toFix );
        unsigned f = isF ? toFix : partner;
        unsigned b = isF ? partner : toFix;

        ++_brokenRelusFixed;

        log( Stringf( "\nAttempting broken-relu fix on var: %s\n", toName( toFix ).ascii() ) );

        double fVal = _assignment[f];
        double bVal = _assignment[b];
        double fDelta;
        double bDelta;

        // f为非负，b为非负
        if ( !FloatUtils::isNegative( fVal ) && !FloatUtils::isNegative( bVal ) )
        {
            fDelta = bVal - fVal;
            bDelta = fVal - bVal;

        } else if ( !FloatUtils::isNegative( fVal ) && FloatUtils::isNegative( bVal )) {
            // f为非负，b为负
            fDelta =  getLeakyValue() * bVal - fVal;
            bDelta = fVal - bVal;

        } else if ( FloatUtils::isNegative( fVal ) && !FloatUtils::isNegative( bVal )) {
            // f为负，b为非负
            fDelta = bVal - fVal;
            bDelta = (1 / getLeakyValue()) * fVal - bVal;

        }else if (FloatUtils::isNegative( fVal ) && FloatUtils::isNegative( bVal )) {
            // f为负，b为负
            fDelta = getLeakyValue() * bVal - fVal;
            bDelta = (1 / getLeakyValue()) * fVal - bVal;

        }else {
            exit( 1 );                           // Unreachable
        }

        bool increaseB = FloatUtils::isPositive( bDelta );
        bool increaseF = FloatUtils::isPositive( fDelta );

        // Always try to fix B first, and if impossible fix f.

        if ( !fixBrokenReluVariable( b, increaseB, bDelta, _brokenReluFixB ) )
            return fixBrokenReluVariable( f, increaseF, fDelta, _brokenReluFixF );

//        if ( !fixBrokenReluVariable( f, increaseF, fDelta, _brokenReluFixF )  )
//            return fixBrokenReluVariable( b, increaseB, bDelta, _brokenReluFixB );

        return true;
    }
    // 先设置好fix时不同variable需要的参数，然后进行fix，优先fix b与f一致，其次再考虑设置f与b一致
    bool fixBrokenRelu_temp( unsigned toFix )
    {
        bool isF = _reluPairs.isF( toFix );
        unsigned partner = _reluPairs.toPartner( toFix );
        unsigned f = isF ? toFix : partner;
        unsigned b = isF ? partner : toFix;

        ++_brokenRelusFixed;

        log( Stringf( "\nAttempting broken-relu fix on var: %s\n", toName( toFix ).ascii() ) );

        double fVal = _assignment[f];
        double bVal = _assignment[b];
        double fDelta;
        double bDelta;

        // 如果forward变量为正，且backward变量为非正，则若要纠正f,则
        if ( FloatUtils::isPositive( fVal ) && !FloatUtils::isPositive( bVal ) )
        {
            fDelta = -fVal;
            bDelta = fVal - bVal;
        }// 如果f为正，且b为正，则各自delta为两者的差值
        else if ( FloatUtils::isPositive( fVal ) && FloatUtils::isPositive( bVal ) )
        {
            fDelta = bVal - fVal;
            bDelta = fVal - bVal;
        }// 如果f为0，b为正数，则纠正f时将f置为b，纠正b时将b置为负数
        else if ( FloatUtils::isZero( fVal ) && FloatUtils::isPositive( bVal ) )
        {
            fDelta = bVal;
            bDelta = -bVal;
        }
        else
            exit( 1 );                           // Unreachable

        bool increaseB = FloatUtils::isPositive( bDelta );
        bool increaseF = FloatUtils::isPositive( fDelta );

        // Always try to fix B first, and if impossible fix f.
        if ( !fixBrokenReluVariable( b, increaseB, bDelta, _brokenReluFixB ) )
            return fixBrokenReluVariable( f, increaseF, fDelta, _brokenReluFixF );

        return true;
    }


    bool fixBrokenReluVariable( unsigned var, bool increase, double &delta, unsigned &_brokenReluStat )
    {
        log( Stringf( "fixBrokenReluVariable Starting: var = %s, delta = %lf\n", toName( var ).ascii(), delta ) );
        //printf( "fixBrokenReluVariable Starting: var = %s, delta = %lf\n", toName( var ).ascii(), delta );

        // 若为non-basic变量
        if ( !_basicVariables.exists( var ) )
        {
            DEBUG(
                    if ( !allVarsWithinBounds() )
                    {
                        printf( "Error! Should not be broken a relu var when we have OOB vars!\n" );
                        exit( 1 );
                    }

                    if ( !canAddToNonBasic( var, delta ) )
                    {
                        // Error: this should not happen. We should not be fixing a broken relu unless
                        // both b and f are within bounds.
                        printf( "Error: var %s is not basic, but can't add delta = %lf to it!\n",
                                toName( var ).ascii(), delta );
                        throw Error( Error::CANT_FIX_BROKEN_RELU, "Unreachable code" );
                    }
            );

            ++_brokenReluStat;
            ++_brokenReluFixByUpdate;

            log( Stringf( "Var %s isn't basic; no pivot needed, simply updating\n", toName( var ).ascii() ) );
            update( var, delta, true );
            return true;
        }// 若为 basic 变量，需要先pivot()，然后再update
        else
        {
            ++_brokenReluStat;
            ++_brokenReluFixByPivot;

            unsigned pivotCandidate;
            if ( !findPivotCandidate( var, increase, pivotCandidate ) ){
//                printf("\n&&&&&&&&&&&& find pivot candidate fail~~~~~~~~~\n");
                return false;
            }

            log( Stringf( "\nPivotAndUpdate: <%s, %5.2lf, %s>\n",
                          toName( var ).ascii(), delta, toName( pivotCandidate ).ascii() ) );

            DEBUG(
                    if ( outOfBounds( var ) )
                    {
                        printf( "Error! Performing a RELU fix when we have an OOB variable\n" );
                        exit( 1 );
                    }
            );

            // Execute the pivot-and-update operation.
            pivot( pivotCandidate, var );
            update( var, delta, true );

            return true;
        }
    }

    void turnAlmostZeroToZero( double &x )
    {
        if ( FloatUtils::isZero( x ) )
            x = 0.0;
    }

    void update( unsigned variable, double delta, bool ignoreRelu = false )
    {
        // delta可能是正数，也可能是负数，

        if ( FloatUtils::isZero( delta ) )
            return;

        log( Stringf( "\t\tUpdate: %s += %.2lf\n", toName( variable ).ascii(), delta ) );

        // 更新值操作
        _assignment[variable] += delta;
        turnAlmostZeroToZero( _assignment[variable] );
        computeVariableStatus( variable );

        const Tableau::Entry *columnEntry = _tableau.getColumn( variable );	// variable is i
        const Tableau::Entry *current;

        // 更新了non-basic变量之后，要根据Tableau里的系数，将凡是有这个non-basic出现的式子的basic的值一起更新，
        // 每次赋值之后，都要对接近0的变量进行置0操作，并根据当前assignment重新计算当前的变量的状态
        while ( columnEntry != NULL )
        {
            current = columnEntry;      // current is the non-basic in this time to update
            columnEntry = columnEntry->nextInColumn();  // columnEntry is the same non-basic which appear in the next row

            unsigned row = current->getRow();   // get the row number ,in order to compare with the variable, if it equals, the value is always -1, no need to update
            if ( row != variable )
            {
                _assignment[row] += delta * current->getValue();
                turnAlmostZeroToZero( _assignment[row] );
                computeVariableStatus( row );
            }
        }

        // If the updated variable was Relu, we might need to fix the relu invariant
        // 如果被更新的是relu变量，那么还要考虑是否需要fix
        // ignoreRelu决定了是否需要忽略relu关系，默认为false,则需要考虑（initialUpdate中没有传值，使用的是默认值，所以还是需要考虑
        if ( _reluPairs.isRelu( variable ) && !ignoreRelu )
        {
            unsigned partner = _reluPairs.toPartner( variable );
            bool variableIsF = _reluPairs.isF( variable );
            unsigned b = variableIsF ? partner : variable;
            unsigned f = variableIsF ? variable : partner;

            log( Stringf( "Update was on relu. Parnter = %u\n", partner ) );

            // If the partner is basic, it's okay for the pair to be broken
            // 现在variable肯定不是basic, 所以如果它的partner是basic , 可以允许暂时broken,
            // because the rules made by algorithm need this, if it is basic , it might be changed by another non-basic
            // only if all the relu pair's two variable are non-baisc , we need to fix them

            if ( _basicVariables.exists( partner ) )
            {
                log( "Partner is basic. ignoring...\n" );
                return;
            }

            // The partner is NOT basic. If the connection is broken,
            // we can fix the partner, if needed.

            // 但如果partner是non-basic, 那么下一步就需要fix这个partner

            log( "Parnter is NOT basic. Checking if more work is needed...\n" );

            log( Stringf( "b = %u, f = %u, bVal = %lf, fVal = %lf\n", b, f, _assignment[b], _assignment[f] ) );

            // _dissolvedReluVariables是Map类型，表示已经被eliminate的ReluPair

            // 如果f存在于_dissolvedReluVariables，那就不需要做任何操作，直接返回
            if ( _dissolvedReluVariables.exists( f ) )
            {
                log( "Pair has been disolved, don't care about a violation\n" );
                return;
            }

            // 如果没有broken，那就不需要任何操作，直接返回
            if ( !reluPairIsBroken( b, f ) )
            {
                log( "relu pair is NOT broken\n" );
                return;
            }

            // 排除了：1、目前已经赋值的variable的partner是基变量
            //		2、forward变量已被eliminate
            //		3、reluPairs没有broken的情况，下面就该进行fix操作

            // 首先判断我们之前assignment的variable是F还是B，如果是F，则partner是B，我们需要fixB,运用updateB规则
            // 否则，partner就是F，运用updateF规则

            if ( variableIsF )
            {
                // We need to fix B. This means setting the value of B to that of F.
                log( Stringf( "Cascading update: fixing non-basic relu partner b = %u\n", b ) );

                // fix backward 变量，

                //由于 forward变量的范围是0~正数，都在backward可取的范围内，所以可以直接赋值
                // 如paper中updateB规则：令B加上其与F的差值（这个差值可能为正，可能为负
                // 调用update函数时，第3个参数传入true,表示这已经是对relu变量进行操作，不需要再进入判断流程

                update( b, _assignment[f] - _assignment[b], true );

                DEBUG(
                        if ( _varToStatus[b] == ABOVE_UB || _varToStatus[b] == BELOW_LB )
                            throw Error( Error::NONBASIC_OUT_OF_BOUNDS,
                                         "After a cascaded b-update, b is non-basic and OOB" );
                );
            }
            else
            {
                // We need to fix F. This means setting the value of F to 0 if B is negative,
                // and otherwise just setting it to B.
                log( Stringf( "Cascading update: fixing non-basic relu partner f = %u\n", f ) );

                // fix forward变量，

                //不能直接赋值，需要根据backward是否大于0，来决定是直接赋值，还是置0

                if ( FloatUtils::isNegative( _assignment[b] ) )
                    update( f, -_assignment[f], true );
                else
                    update( f, _assignment[b] - _assignment[f], true );

                DEBUG(
                        if ( _varToStatus[f] == ABOVE_UB || _varToStatus[f] == BELOW_LB )
                            throw Error( Error::NONBASIC_OUT_OF_BOUNDS,
                                         "After a cascaded f-update, f is non-basic and OOB" );
                );
            }
        }
    }

    void pivot( unsigned nonBasic, unsigned basic )
    {
        ++_numPivots;

        log( Stringf( "\t\tPivot: %s <--> %s\n", toName( basic ).ascii(), toName( nonBasic ).ascii() ) );

        // Sanity checks:
        //检查传入的non-basic和basic是否是正确的
        if ( _basicVariables.exists( nonBasic ) )
            throw Error( Error::ILLEGAL_PIVOT_OP,
                         Stringf( "Non-basic variable %s is basic", toName( nonBasic ).ascii() ).ascii() );

        if ( !_basicVariables.exists( basic ) )
            throw Error( Error::ILLEGAL_PIVOT_OP, "Basic variable isn't basic" );

        // 1、更改_basicVariables里存储的值
        _basicVariables.erase( basic );
        _basicVariables.insert( nonBasic );

        timeval start = Time::sampleMicro();
        unsigned numCalcs = 0;

        // 取得tableau中对应表示两个变量的关系的系数，如果没有，则返回0.0
        double cell = _tableau.getCell( basic, nonBasic );
        double absWeight = FloatUtils::abs( cell );

        if ( FloatUtils::lt( absWeight, NUMBERICAL_INSTABILITY_CONSTANT ) )
        {
            printf( "--- Numerical Instability Warning!! Weight = %.15lf ---\n", absWeight );
        }

        // 将新basic写入Tableau
        _tableau.addScaledRow( basic,
                               ( -1.0 ) / cell,
                               nonBasic,
                // Guarantee a -1 in the (nb,nb) cell
                               nonBasic, -1.0,
                //
                               &numCalcs
        );
        // 将原basic从Tableau中删除，主要是修正该Entry出现在的列链表，删除Enrty,最后删除_rows数组中存储的对应表头，即删除该双向链表
        _tableau.eraseRow( basic );

        log( Stringf( "\t\t\tPivot--clearing %u column entries--starting\n",
                      _tableau.getColumnSize( nonBasic ) ) );

        // 完成当前等式的pivot之后，还要将新Basic出现过的其他basic等式也一起做更改
        const Tableau::Entry *columnEntry = _tableau.getColumn( nonBasic );
        const Tableau::Entry *current;

        while ( columnEntry != NULL )
        {
            current = columnEntry;
            columnEntry = columnEntry->nextInColumn();

            // 如果在删除了原basic行之后，新basic所在的列还有其他值，也就是它还出现在了其他等式的右边，此时必须要进行相关处理，将它所在列的其他行置0，
            // 同时处理置0所在行的其他系数 (通过addScaledRow一并完成）
            if ( current->getRow() != nonBasic )
            {
                _tableau.addScaledRow( nonBasic,
                                       current->getValue(),
                                       current->getRow(),
                        // Guarantee a 0 in the (*,nb) cells
                                       nonBasic, 0.0,
                        //
                                       &numCalcs
                );
            }
        }

        timeval end = Time::sampleMicro();
        log( Stringf( "\t\t\tPivot--clearing column entries--done (Pivot: %u milli, %u calcs)\n",
                      Time::timePassed( start, end ), numCalcs ) );

        _totalPivotTimeMilli += Time::timePassed( start, end );
        _totalPivotCalculationCount += numCalcs;
    }

    void dump()
    {
        if ( !_dumpStates )
            return;

        log( "\nVisiting state:\n" );

        log( "\n" );
        log( "       | " );
        for ( unsigned i = 0; i < _numVariables; ++i )
            log( Stringf( "%6s", toName( i ).ascii() ) );
        log( " | Assignment               " );
        log( "\n" );
        for ( unsigned i = 0; i < 9 + ( _numVariables * 6 ) + 13 + 15; ++i )
            log( "-" );
        log( "\n" );

        for ( unsigned i = 0; i < _numVariables; ++i )
        {
            if ( _basicVariables.exists( i ) )
                log( " B " );
            else
                log( "   " );

            log( Stringf( "%4s| ", toName( i ).ascii() ) );
            for ( unsigned j = 0; j < _numVariables; ++j )
            {
                if ( !FloatUtils::isZero( _tableau.getCell( i, j ) ) )
                    log( Stringf( "%6.2lf", _tableau.getCell( i, j ) ) );
                else
                    log( "      " );
            }
            log( " | " );

            if ( _lowerBounds[i].finite() )
                log( Stringf( "%5.2lf <= ", _lowerBounds[i].getBound() ) );
            else
                log( "         " );

            log( Stringf( "%.4lf", _assignment[i] ) );

            if ( outOfBounds( i ) || ( activeReluVariable( i ) && partOfBrokenRelu( i ) ) )
                log( " * " );
            else
                log( "   " );

            if ( _upperBounds[i].finite() )
                log( Stringf( "<= %5.2lf", _upperBounds[i].getBound() ) );
            else
                log( "         " );

            log( "\n" );
        }

        log( "\n" );
    }

    bool canAddToNonBasic( unsigned variable, double delta )
    {
        if ( FloatUtils::isZero( delta ) )
            return true;

        bool positive = FloatUtils::isPositive( delta );
        VariableStatus status = _varToStatus[variable];

        if ( status == ABOVE_UB || status == BELOW_LB )
            throw Error( Error::NONBASIC_OUT_OF_BOUNDS );

        if ( status == FIXED )
            return false;

        if ( positive )
        {
            if ( status == AT_UB && FloatUtils::gt( delta, OOB_EPSILON ) )
                return false;

            if ( !_upperBounds[variable].finite() )
                return true;

            return FloatUtils::lte( _assignment[variable] + delta, _upperBounds[variable].getBound(), OOB_EPSILON );
        }
        else
        {
            if ( status == AT_LB && FloatUtils::lt( delta, -OOB_EPSILON ) )
                return false;

            if ( !_lowerBounds[variable].finite() )
                return true;

            return FloatUtils::gte( _assignment[variable] + delta, _lowerBounds[variable].getBound(), OOB_EPSILON );
        }
    }

    bool tooLow( unsigned variable ) const
    {
        return _varToStatus.at( variable ) == VariableStatus::BELOW_LB;
    }

    bool canDecrease( unsigned variable ) const
    {
        return
                _varToStatus.at( variable ) == VariableStatus::BETWEEN ||
                _varToStatus.at( variable ) == VariableStatus::AT_UB ||
                _varToStatus.at( variable ) == VariableStatus::ABOVE_UB;
    }

    bool tooHigh( unsigned variable ) const
    {
        return _varToStatus.at( variable ) == VariableStatus::ABOVE_UB;
    }

    bool canIncrease( unsigned variable ) const
    {
        return
                _varToStatus.at( variable ) == VariableStatus::BETWEEN ||
                _varToStatus.at( variable ) == VariableStatus::AT_LB ||
                _varToStatus.at( variable ) == VariableStatus::BELOW_LB;
    }

    bool outOfBounds( unsigned variable ) const
    {
        return tooLow( variable ) || tooHigh( variable );
    }

    void log( String message )
    {
        if ( _logging )
            printf( "%s", message.ascii() );
    }

    double getAssignment( unsigned variable ) const
    {
        return _assignment[variable];
    }

    void setLogging( bool value )
    {
        _logging = value;
    }

    void setDumpStates( bool value )
    {
        _dumpStates = value;
    }

    unsigned numStatesExplored() const
    {
        return _numCallsToProgress;
    }

    bool fixedAtZero( unsigned var ) const
    {
        return
                ( _varToStatus.get( var ) == VariableStatus::FIXED ) &&
                FloatUtils::isZero( _upperBounds[var].getBound() );
    }

    bool eliminateAuxVariables()
    {
        log( "eliminateAuxVariables starting\n" );
        computeVariableStatus();    // 设置_varToStatus

        Set<unsigned> initialAuxVariables = _basicVariables;

        for ( const auto &aux : initialAuxVariables )
        {
            // eliminateIfPossible只有在没找到可供pivot的候选者时才会返回false,
            // 否则只要有过pivot,不论是否eliminate，都会返回true
            // 只要有其中一个辅助变量没法再pivot,就返回false
            if ( !eliminateIfPossible( aux ) )
            {
                log( "eliminateAuxVariables finished UNsuccessfully\n" );
                return false;
            }
        }

        log( "eliminateAuxVariables finished successfully\n" );
        return true;
    }

    bool eliminateIfPossible( unsigned var )
    {
        //初始时，basic都是另外引入的辅助变量，所以这里只可能是辅助变量，不可能是relu，
        // 如果此函数只在初始时调用，则没问题，如果后续运行中继续调用，可能会有问题
        if ( _reluPairs.isRelu( var ) )
        {
            printf( "Attempted to eliminate a relu variable. They shouldn't be marked as aux\n" );
            exit( 1 );
        }
        // 判断是否超越下界,这里的var只会是最初调用时的辅助变量，所以一定时候basic
        bool increase = tooLow( var );

        // 根据越界情况，取delta为上或下界与value的差值
        double delta = increase?
                       ( _lowerBounds[var].getBound() - _assignment[var] ) :
                       ( _upperBounds[var].getBound() - _assignment[var] );

        // 查找可供进行pivot的候选者
        // 没有找到，返回false，退出函数
        unsigned pivotCandidate;
        if ( !findPivotCandidate( var, increase, pivotCandidate, false ) )
            return false;

        log( Stringf( "\nPivotAndUpdate: <%s, %5.2lf, %s>\n",
                      toName( var ).ascii(),
                      delta,
                      toName( pivotCandidate ).ascii() ) );

        // 进行pivot操作
        pivot( pivotCandidate, var );
        // 更新,如果delta为0，则自动返回
        update( var, delta );

        if ( !fixedAtZero( var ) )
        {
            printf( "eliminateIfPossible called for a non fixed-at-zero variable\n" );
            return true;
        }

        log( Stringf( "\nVariable %s fixed at zero. Eliminating...\n", toName( var ).ascii() ) );
        _tableau.eraseColumn( var );
        _eliminatedVars.insert( var );
        ++_numEliminatedVars;

        return true;
    }

    bool findPivotCandidate( unsigned variable, bool increase, unsigned &pivotCandidate,
                             bool ensureNumericalStability = true )
    {
        const Tableau::Entry *rowEntry = _tableau.getRow( variable );
        const Tableau::Entry *current;

        unsigned column;

        bool found = false;
        unsigned leastEvilNonBasic = 0;
        double leastEvilWeight = 0.0;

        while ( rowEntry != NULL )
        {
            current = rowEntry;
            rowEntry = rowEntry->nextInRow();

            column = current->getColumn();

            // Ignore self
            if ( column == variable )
                continue;


            const double coefficient = current->getValue();
            bool positive = FloatUtils::isPositive( coefficient );

            // increase表示是否越下界，若是则为true ,需要增加 , basic变量BELOW_LB
            // !increase表示不越下界，那么就是 AT_LB,BETWEEN，AT_UB,ABOVE_UB
            // canIncrease()表示一定是小于上界的，不会大于也不会等于，根据变量此时的状态判断是否可以再继续增大值，如果是越下界BELOW_LB、等于下界AT_LB、上下界之间BETWEEN，则返回ture
            // canDecrease()表示一定是大于下界的，不会大于也不会等于，与上相反，如果是越上界ABOVE_UB，等于上界AT_UB、上下界之间BETWEEN，则返回true

            // 总结：1、原basic越下界<、系数为正、且这个变量的value可以被增加
            // 2、原basic越下界<、系数为负，且这个变量的value可以被减少
            // 3、原basic非越下界>=（大于等于下界都可以，可能越上界，也可能不越）、系数为正、value可以被减少
            // 4、原basic非越下界>=、系数为负，value可以被增加
            // 有以上4种情况的变量，可以进入下一步，否则就意味着当前变量不满足paper中slack的要求，退出此次循环查找下一个

            // 这里其实只要大致方向对就OK，即以下界为分割点，只要小于下界，那么就只能进行加正和减负操作
            // 只要大于等于下界，那么就可以进行减正和加负操作，
            // 我们并不保证操作后的non-basic(即将变成新basic)一定在范围内，
            // 如果小于下界，可能加正、减负之后还是小于下界，如果大于等于下界，可能减正、加负之后会小于下界或仍然大于上界，
            // 但此时并不考虑这些，只考虑运算的演进方向是对的
            // 所以如果初始化时，一个辅助变量在左边的等式，如果找不到可以变换的值，那么单纯形法就要报错，因为最终辅助变量都是要变换到右边设值为0的

            if ( !( ( increase && ( positive ) && canIncrease( column ) ) ||
                    ( increase && ( !positive ) && canDecrease( column ) ) ||
                    ( !increase && ( positive ) && canDecrease( column ) ) ||
                    ( !increase && ( !positive ) && canIncrease( column ) ) ) )
            {
                // The variable does not fit the direction we need.
                continue;
            }

            double weight = FloatUtils::abs( coefficient );

            // 默认值为true,但是传入值为false,一定会进入if,返回true
            // ensureNumericalStability是指是否需要保证数字的稳定性，因为当数字太小时，由于计算机固有误差，会近似等于0
            // 而此时刚刚开始进行计算，传入false,可以忽略这一要求，保证在满足slack条件的情况下一定能找到pivot候选者
            if ( !ensureNumericalStability || FloatUtils::gte( weight, NUMBERICAL_INSTABILITY_CONSTANT ) )
            {
                pivotCandidate = column;
                return true;
            }

            // Have a candidate with a small pivot coefficient
            // 如果找到了一个候选者，但是它的值非常小，先将第一个记录下来，再进行后续比较，如果后续还发现了符合条件，而权重值更大的pivot候选者，就更新least记录
            found = true;
            if ( FloatUtils::gt( weight, leastEvilWeight ) )
            {
                leastEvilWeight = weight;
                leastEvilNonBasic = column;
            }
        }
        // 在遍历完所有变量之后，least中记录的是权重值最大的候选者，将其返回
        if ( found )
        {
            log( Stringf( "findPivotCandidate: forced to pick a bad candidate! Weight = %lf\n", leastEvilWeight ) );
            pivotCandidate = leastEvilNonBasic;
            return true;
        }

        // 没找到pivot候选者
        return false;
    }

    const VariableBound *getLowerBounds() const
    {
        return _lowerBounds;
    }

    const VariableBound *getUpperBounds() const
    {
        return _upperBounds;
    }

    double getLowerBound( unsigned var ) const
    {
        DEBUG(
                if ( !_lowerBounds[var].finite() )
                    throw Error( Error::LOWER_BOUND_IS_INFINITE );
        );

        return _lowerBounds[var].getBound();
    }

    double getUpperBound( unsigned var ) const
    {
        DEBUG(
                if ( !_upperBounds[var].finite() )
                    throw Error( Error::UPPER_BOUND_IS_INFINITE );
        );

        return _upperBounds[var].getBound();
    }

    Set<unsigned> getBasicVariables() const
    {
        return _basicVariables;
    }

    ReluPairs *getReluPairs()
    {
        return &_reluPairs;
    }

    void setLowerBounds( const List<VariableBound> &lowerBounds )
    {
        unsigned i = 0;
        for ( const auto &it : lowerBounds )
        {
            _lowerBounds[i] = it;
            ++i;
        }
    }

    void setUpperBounds( const List<VariableBound> &upperBounds )
    {
        unsigned i = 0;
        for ( const auto &it : upperBounds )
        {
            _upperBounds[i] = it;
            ++i;
        }
    }

    void setBasicVariables( const Set<unsigned> &basicVariables )
    {
        _basicVariables = basicVariables;
    }

    void setReluPairs( const ReluPairs &reluPairs )
    {
        _reluPairs = reluPairs;
    }

    void backupIntoMatrix( Tableau *matrix ) const
    {
        _tableau.backupIntoMatrix( matrix );
    }

    void restoreFromMatrix( Tableau *matrix )
    {
        matrix->backupIntoMatrix( &_tableau );

        DEBUG(
                log( "Printing matrix after restoration\n" );
                log( "****\n" );
                dump();
                log( "****\n\n" );
        );
    }

    const double *getAssignment() const
    {
        return _assignment;
    }

    void setAssignment( const List<double> &assignment )
    {
        unsigned i = 0;
        for ( const auto &it : assignment )
        {
            _assignment[i] = it;
            ++i;
        }
    }

    void makeAllBoundsFinite()
    {
        // 计算上或下界中有无穷大值的变量的个数
        countVarsWithInfiniteBounds();
        log( Stringf( "makeAllBoundsFinite -- Starting (%u vars with infinite bounds)\n", _varsWithInfiniteBounds ) );
//        printf("\n----- printStatistics() ----\n");
//        printStatistics();

        // 不懂
        dump();
        for ( const auto &basic : _basicVariables )
            makeAllBoundsFiniteOnRow( basic );

        countVarsWithInfiniteBounds();
        log( Stringf( "makeAllBoundsFinite -- Done (%u vars with infinite bounds)\n", _varsWithInfiniteBounds ) );
        printBounds();

//        printf("\n----- printStatistics() ----\n");
//        printStatistics();

        if ( _varsWithInfiniteBounds != 0 )
            throw Error( Error::EXPECTED_NO_INFINITE_VARS );
    }

    void makeAllBoundsFiniteOnRow( unsigned basic )
    {
//        printf("~~~~basic: %u\n", basic);
        // 取得当前basic变量的row入口
        const Tableau::Entry *row = _tableau.getRow( basic );
        const Tableau::Entry *tighteningVar = NULL;

        // 遍历查找具有无穷大界限的那个变量，地址存储在tighteningVar中，
        // 如果有多个，报错
        while ( row != NULL )
        {
//            printf("~~~~row->getColumn(): %u\n", row->getColumn());

            if ( !_upperBounds[row->getColumn()].finite() || !_lowerBounds[row->getColumn()].finite() )
            {
                if ( tighteningVar != NULL )
                    throw Error( Error::MULTIPLE_INFINITE_VARS_ON_ROW );

                tighteningVar = row;
            }

            row = row->nextInRow();
        }

        // It's possible that there are no infinite vars on this row - e.g., if the user supplied
        // bounds on the outputs.
        if ( !tighteningVar )
            return;

        unsigned tighteningVarIndex = tighteningVar->getColumn();

        // 如果找到了无穷大界限的变量，设置缩放量为 -1.0 / 888，
        double scale = -1.0 / tighteningVar->getValue();

        row = _tableau.getRow( basic );
        const Tableau::Entry *current;

        double max = 0.0;
        double min = 0.0;

        // 同上方法再次遍历，
        while ( row != NULL )
        {
            current = row;
            row = row->nextInRow();

            if ( current->getColumn() == tighteningVarIndex )
                continue;

            double coefficient = current->getValue() * scale;
            if ( FloatUtils::isPositive( coefficient ) )
            {
                max += _upperBounds[current->getColumn()].getBound() * coefficient;
                min += _lowerBounds[current->getColumn()].getBound() * coefficient;
            }
            else
            {
                min += _upperBounds[current->getColumn()].getBound() * coefficient;
                max += _lowerBounds[current->getColumn()].getBound() * coefficient;
            }
        }

        // 如果是上界无穷大，或者计算出的max小于上界，即可以缩小界限的上界范围，则更新上界
        if ( !_upperBounds[tighteningVarIndex].finite() ||
             FloatUtils::lt( max, _upperBounds[tighteningVarIndex].getBound() ) )
            updateUpperBound( tighteningVarIndex, max, 0 );

        //  如果是下界无穷大，或者计算出的min大于下界，即可以缩小下界的范围。
        if ( !_lowerBounds[tighteningVarIndex].finite() ||
             FloatUtils::gt( min, _lowerBounds[tighteningVarIndex].getBound() ) )
            updateLowerBound( tighteningVarIndex, min, 0 );

        computeVariableStatus( tighteningVarIndex );
        if ( !_basicVariables.exists( tighteningVarIndex ) && outOfBounds( tighteningVarIndex ) )
            update( tighteningVarIndex,
                    _lowerBounds[tighteningVarIndex].getBound() - _assignment[tighteningVarIndex] );
    }

    void setUseApproximation( bool value )
    {
        _useApproximations = value;
    }

    void setFindAllPivotCandidates( bool value )
    {
        _findAllPivotCandidates = value;
    }

    void showDissolvedMergeReluPairs(){
//        printf("~~~~~When finished, the relu variable be merged are the following:\n");
        for (unsigned i = 0; i < _numVariables; i++) {
            if (isDissolvedBVariable(i)) {
//                printf("~~~~~ %s : %u\n ", toName( i ).ascii(), i);
            }
        }
    }
    // 判断某个变量是否是已经被求解了的，
    bool isDissolvedBVariable( unsigned variable ) const
    {
        // 首先需要是一个relu变量
        if ( !_reluPairs.isRelu( variable ) )
            return false;

        // 其次要是个b变量
        if ( _reluPairs.isF( variable ) )
            return false;

        unsigned f = _reluPairs.toPartner( variable );

        // 取得b变量的f，用来判断是否已存在于_dissolvedReluVariables里面，如果不存在，那么返回false表示没有被求解
        if ( !_dissolvedReluVariables.exists( f ) )
            return false;

        // 如果确实已经被添加到_dissolvedReluVariables里面，那么看它的类型是不是merge，如果是,返回true
        return _dissolvedReluVariables.at( f ) == TYPE_MERGE;
    }

    bool isEliminatedVar( unsigned variable ) const
    {
        return _eliminatedVars.exists( variable );
    }

    void markReluVariableDissolved( unsigned variable, ReluDissolutionType type )
    {
        log( Stringf( "Mark var as dissolved: %u (Type: %s)\n",
                      variable,
                      type == TYPE_SPLIT ? "Split" : "Merge" ) );

        DEBUG(
                if ( _dissolvedReluVariables.exists( variable ) )
                {
                    printf( "Error -- this variable was already marked as dissolved!\n" );
                    exit( 1 );
                }
        );

        _dissolvedReluVariables[variable] = type;
    }

    void incNumSplits()
    {
        ++_numStackSplits;
    }

    void incNumMerges()
    {
        ++_numStackMerges;
    }

    void incNumPops()
    {
        ++_numStackPops;
    }

    void incNumStackVisitedStates()
    {
        ++_numStackVisitedStates;
    }

    void setCurrentStackDepth( unsigned depth )
    {
        _currentStackDepth = depth;
        if ( _currentStackDepth > _maximalStackDepth )
            _maximalStackDepth = _currentStackDepth;
    }

    void setMinStackSecondPhase( unsigned depth )
    {
        if ( ( depth < _minStackSecondPhase ) || ( _minStackSecondPhase == 0 ) )
            _minStackSecondPhase = depth;
    }

    unsigned getColumnSize( unsigned column ) const
    {
        return _tableau.getColumnSize( column );
    }

    Map<unsigned, ReluDissolutionType> getDissolvedReluPairs() const
    {
        return _dissolvedReluVariables;
    }

    void setDissolvedReluPairs( const Map<unsigned, ReluDissolutionType> &pairs )
    {
        _dissolvedReluVariables = pairs;
    }

    unsigned reluVarToF( unsigned variable ) const
    {
        if ( _reluPairs.isF( variable ) )
            return variable;

        return _reluPairs.toPartner( variable );
    }

    bool isReluVariable( unsigned variable ) const
    {
        return _reluPairs.isRelu( variable );
    }

    void printAssignment()
    {
        printf( "\nCurrent assignment:\n" );
        for ( unsigned i = 0; i < _numVariables; ++i )
        {
            if ( _eliminatedVars.exists( i ) || !outOfBounds( i ) )
                continue;

            printf( "\t%u: %.10lf <= %.10lf <= %.10lf",
                    i, _lowerBounds[i].getBound(), _assignment[i], _upperBounds[i].getBound() );
            if ( outOfBounds( i ) )
                printf( "  ***" );
            if ( _basicVariables.exists( i ) )
                printf( " B" );
            printf( "\n" );
        }
        printf( "\n" );
    }

    unsigned getNumVariables() const
    {
        return _numVariables;
    }

    const Tableau::Entry *getColumn( unsigned column ) const
    {
        return _tableau.getColumn( column );
    }

    const Tableau::Entry *getRow( unsigned row ) const
    {
        return _tableau.getRow( row );
    }

    double getCell( unsigned row, unsigned column ) const
    {
        return _tableau.getCell( row, column );
    }

    Set<unsigned> getEliminatedVars() const
    {
        return _eliminatedVars;
    }

    const Tableau *getTableau() const
    {
        return &_tableau;
    }

    void storePreprocessedMatrix()
    {
        checkInvariants();

        // 将当前初始变换后的Tableau存入_preprocessedTableau
        // Tableau //Reluplex构造方法中初始化，initializeCell中添加，pivot中修改
        _tableau.backupIntoMatrix( &_preprocessedTableau );

        // 初始化时暂时不涉及
        _preprocessedDissolvedRelus = _dissolvedReluVariables;


        // _basicVariables;	//markBasic中初始化，pivot中修改
        _preprocessedBasicVariables = _basicVariables;

        // 存储当前的赋值到 _preprocessedAssignment
        // _assignment	// Reluplex构造方法中默认初始化为0，update中修改
        memcpy( _preprocessedAssignment, _assignment, sizeof(double) * _numVariables );

        // // setUpperBound中初始化，makeAllBoundsFiniteOnRow > updateUpperBound中修改
        for ( unsigned i = 0; i < _numVariables; ++i )
        {
            // It is assumed that these bounds are all at level 0.
            _preprocessedLowerBounds[i] = _lowerBounds[i];
            _preprocessedUpperBounds[i] = _upperBounds[i];
        }
    }

    void restoreTableauFromBackup( bool keepCurrentBasicVariables = true )
    {

        printf("\n@@@@@@@@@@@@ restore from restoreTableauFromBackup: _preprocessedDissolvedRelus.size : %u\n", _preprocessedDissolvedRelus.size());

        timeval start = Time::sampleMicro();

        ++_numberOfRestorations;

        printf( "\n\n\t\t !!! Restore tableau from backup starting !!!\n" );
        double *backupLowerBounds = new double[_numVariables];
        double *backupUpperBounds = new double[_numVariables];
        unsigned *backupLowerBoundLevels = new unsigned[_numVariables];
        unsigned *backupUpperBoundLevels = new unsigned[_numVariables];

        Set<unsigned> backupBasicVariables = _basicVariables;

        for ( unsigned i = 0; i < _numVariables; ++i )
        {
            DEBUG({
                      if ( FloatUtils::lt( _lowerBounds[i].getBound(), _preprocessedLowerBounds[i].getBound() ) )
                      {
                          printf( "Error with a decreasing LB\n" );
                          exit( 1 );
                      }

                      if ( FloatUtils::gt( _upperBounds[i].getBound(), _preprocessedUpperBounds[i].getBound() ) )
                      {
                          printf( "Error with an increasing UBs\n" );
                          exit( 1 );
                      }

                      if ( FloatUtils::gt( _lowerBounds[i].getBound(), _upperBounds[i].getBound() ) )
                      {
                          printf( "Error! LB > UB\n" );
                          exit( 1 );
                      }
                  });

            backupLowerBounds[i] = _lowerBounds[i].getBound();
            backupUpperBounds[i] = _upperBounds[i].getBound();
            backupLowerBoundLevels[i] = _lowerBounds[i].getLevel();
            backupUpperBoundLevels[i] = _upperBounds[i].getLevel();

            _lowerBounds[i] = _preprocessedLowerBounds[i];
            _upperBounds[i] = _preprocessedUpperBounds[i];
            _lowerBounds[i].setLevel( 0 );
            _upperBounds[i].setLevel( 0 );
        }

        Map<unsigned, ReluDissolutionType> backupDissolved = _dissolvedReluVariables;

        // First restore the original tableau and the state of relus
        _preprocessedTableau.backupIntoMatrix( &_tableau );
        _dissolvedReluVariables = _preprocessedDissolvedRelus;
        memcpy( _assignment, _preprocessedAssignment, sizeof(double) * _numVariables );
        _basicVariables = _preprocessedBasicVariables;
        computeVariableStatus();

        // Since we restored the original tableau, all invariants should hold.
        checkInvariants();

        // Next, go over the current bounds and assert them, one by one
        for ( unsigned i = 0; i < _numVariables; ++i )
        {
            double newLb = backupLowerBounds[i];
            double newUb = backupUpperBounds[i];

            // If the variable is not an active relu (in the original tableau), just update it.
            if ( !activeReluVariable( i ) )
            {
                if ( ( !_lowerBounds[i].finite() ) || FloatUtils::gt( newLb, _lowerBounds[i].getBound() ) )
                    updateLowerBound( i, newLb, backupLowerBoundLevels[i] );

                if ( ( !_upperBounds[i].finite() ) || FloatUtils::lt( newUb, _upperBounds[i].getBound() ) )
                    updateUpperBound( i, newUb, backupUpperBoundLevels[i] );

                continue;
            }

            // Dealing with active relu variables. Only handle F's.
            if ( !_reluPairs.isF( i ) )
            {
                continue;
            }

            unsigned f = i;
            unsigned b = _reluPairs.toPartner( i );

            double bLower = backupLowerBounds[b];
            double bUpper = backupUpperBounds[b];

            // Now, it matters whether F and B are still active in the tableau we're creating.
            if ( !backupDissolved.exists( f ) )
            {
                // Stil Active. Upper bounds should match and be equal, so can just update one of them and cascade.
                if ( ( !_upperBounds[i].finite() ) || FloatUtils::lt( newUb, _upperBounds[i].getBound() ) )
                    updateUpperBound( f, newUb, backupUpperBoundLevels[f] );

                // The lower bound for F was not updated, or the pair would be merged. Hence, just update B.
                // The update function knows how to handle this propertly.
                if ( ( !_lowerBounds[b].finite() ) || FloatUtils::gt( bLower, _lowerBounds[b].getBound() ) )
                    updateLowerBound( b, bLower, backupLowerBoundLevels[b] );
            }

            else if ( backupDissolved[f] == TYPE_SPLIT )
            {
                // Not active, due to a split. B's upper bound must be non-positive.
                // Update, and this will cause a split.
                if ( ( !_upperBounds[b].finite() ) || FloatUtils::lt( bUpper, _upperBounds[b].getBound() ) )
                    updateUpperBound( b, bUpper, backupUpperBoundLevels[b] );

                // Maybe the pair was broken at an earlier update, so fix F's levels individually.
                _upperBounds[f].setLevel( backupUpperBoundLevels[f] );

                // Also, normally update b's lower bound
                if ( ( !_lowerBounds[b].finite() ) || FloatUtils::gt( bLower, _lowerBounds[b].getBound() ) )
                    updateLowerBound( b, bLower, backupLowerBoundLevels[b] );
            }

            else
            {
                // Not active, due to a merge. B's lower bound must be non-negative.
                // Update, and this will cause the merge.
                if ( ( !_lowerBounds[b].finite() ) || FloatUtils::gt( bLower, _lowerBounds[b].getBound() ) )
                    updateLowerBound( b, bLower, backupLowerBoundLevels[b] );

                // Now that this is done, update both bounds for F - maybe they were tightend later.
                // (I.e., maybe the new lower bound for f is tighter than the one for b).
                if ( ( !_lowerBounds[f].finite() ) || FloatUtils::gt( newLb, _lowerBounds[f].getBound() ) )
                    updateLowerBound( f, newLb, backupLowerBoundLevels[f] );

                if ( ( !_upperBounds[f].finite() ) || FloatUtils::lt( newUb, _upperBounds[f].getBound() ) )
                    updateUpperBound( f, newUb, backupUpperBoundLevels[f] );
            }
        }

        DEBUG({
                  if ( backupDissolved != _dissolvedReluVariables )
                  {
                      printf( "Error - didnt get the same set of dissolved relus\n" );
                      exit( 1 );
                  }

                  checkInvariants();

                  // The current bounds should be the same as before, except for merged B's.
                  for ( unsigned i = 0; i < _numVariables; ++i )
                  {
                      // Ignore merged B's
                      bool isB = _reluPairs.isB( i );
                      if ( isB )
                      {
                          unsigned f = _reluPairs.toPartner( i );
                          if ( _dissolvedReluVariables.exists( f ) &&
                               _dissolvedReluVariables[f] == TYPE_MERGE )
                          {
                              continue;
                          }
                      }

                      // Otherwise, bounds need to be equal.
                      if ( FloatUtils::areDisequal( _lowerBounds[i].getBound(), backupLowerBounds[i] ) )
                      {
                          // This is only allowed if we have a relu pair in which b is eliminated;
                          // its bound may be different than expected

                          if ( _reluPairs.isRelu( i ) )
                          {
                              unsigned f = _reluPairs.isF( i ) ? i : _reluPairs.toPartner( i );
                              unsigned b = _reluPairs.toPartner( f );

                              if ( !( ( i == b ) && ( _tableau.getColumnSize( b ) == 0 ) ) )
                              {
                                  printf( "Error in lower bounds for var %u. %lf != %lf\n",
                                          i, _lowerBounds[i].getBound(), backupLowerBounds[i] );

                                  printf( "Checking. b = %u, f = %u. Dissolved? %s.\n",
                                          b, f, _dissolvedReluVariables.exists( f ) ? "YES" : "NO" );
                                  printf( "B's column size: %u. F's column size: %u\n",
                                          _tableau.getColumnSize( b ), _tableau.getColumnSize( f ) );

                                  printf( "Original bounds for b = %u: lower = %.15lf, upper = %.15lf",
                                          b, backupLowerBounds[b], backupUpperBounds[b] );
                                  printf( "Original bounds for f = %u: lower = %.15lf, upper = %.15lf",
                                          f, backupLowerBounds[f], backupUpperBounds[f] );

                                  printf( "And, bounds after the update:\n" );
                                  printf( "\tb = %u: lower = %.15lf, upper = %.15lf",
                                          b, _lowerBounds[b].getBound(), _upperBounds[b].getBound() );
                                  printf( "\tf = %u: lower = %.15lf, upper = %.15lf",
                                          f, _lowerBounds[f].getBound(), _upperBounds[f].getBound() );

                                  printf( "Not the case of an eliminated b variable!\n" );
                                  exit( 1 );
                              }
                          }
                          else
                          {
                              printf( "Error in lower bounds for var %u. %lf != %lf\n", i, _lowerBounds[i].getBound(), backupLowerBounds[i] );
                              printf( "Not relu!\n" );
                              exit( 1 );
                          }
                      }

                      if ( FloatUtils::areDisequal( _upperBounds[i].getBound(), backupUpperBounds[i] ) )
                      {
                          if ( _reluPairs.isRelu( i ) )
                          {
                              unsigned f = _reluPairs.isF( i ) ? i : _reluPairs.toPartner( i );
                              unsigned b = _reluPairs.toPartner( f );

                              if ( !( ( i == b ) && ( _tableau.getColumnSize( b ) == 0 ) ) )
                              {
                                  printf( "Error in upper bounds for var %u. %.15lf != %.15lf\n",
                                          i, _upperBounds[i].getBound(), backupUpperBounds[i] );

                                  printf( "Checking. b = %u, f = %u. Dissolved? %s.\n",
                                          b, f, _dissolvedReluVariables.exists( f ) ? "YES" : "NO" );
                                  printf( "B's column size: %u. F's column size: %u\n",
                                          _tableau.getColumnSize( b ), _tableau.getColumnSize( f ) );

                                  printf( "Original bounds for b = %u: lower = %.15lf, upper = %.15lf\n",
                                          b, backupLowerBounds[b], backupUpperBounds[b] );
                                  printf( "Original bounds for f = %u: lower = %.15lf, upper = %.15lf\n",
                                          f, backupLowerBounds[f], backupUpperBounds[f] );

                                  printf( "And, bounds after the update:\n" );
                                  printf( "\tb = %u: lower = %.15lf, upper = %.15lf\n",
                                          b, _lowerBounds[b].getBound(), _upperBounds[b].getBound() );
                                  printf( "\tf = %u: lower = %.15lf, upper = %.15lf\n",
                                          f, _lowerBounds[f].getBound(), _upperBounds[f].getBound() );

                                  printf( "Not the case of an eliminated b variable!\n" );
                                  exit( 1 );
                              }
                          }
                          else
                          {
                              printf( "Error in upper bounds for var %u. %.15lf != %.15lf\n",
                                      i, _upperBounds[i].getBound(), backupUpperBounds[i] );
                              printf( "Not relu!\n" );
                              exit( 1 );
                          }
                      }

                      // Levels should also be equal
                      if ( _lowerBounds[i].getLevel() != backupLowerBoundLevels[i] )
                      {
                          printf( "Error restoring lower bound for variable %s. "
                                  "Expected: %u, got: %u\n", toName( i ).ascii(),
                                  backupLowerBoundLevels[i], _lowerBounds[i].getLevel() );
                          exit( 1 );
                      }

                      if ( _upperBounds[i].getLevel() != backupUpperBoundLevels[i] )
                      {
                          printf( "Error restoring upper bound for variable %s. "
                                  "Expected: %u, got: %u\n", toName( i ).ascii(),
                                  backupLowerBoundLevels[i], _lowerBounds[i].getLevel() );
                          exit( 1 );
                      }
                  }
              });

        // The tableau now has the proper bounds and proper eliminated variables.
        // All that remains is to pivot to the correct basis.
        if ( keepCurrentBasicVariables )
        {
            printf( "\t\t\tRestoring basics\n" );
            Set<unsigned> shouldBeBasic = Set<unsigned>::difference( backupBasicVariables, _basicVariables );
            Set<unsigned> shouldntBeBasic = Set<unsigned>::difference( _basicVariables, backupBasicVariables );

            adjustBasicVariables( shouldBeBasic, shouldntBeBasic );
        }
        else
        {
            printf( "\t\t\tNot restoring basics\n" );
        }

        DEBUG( checkInvariants(); );

        delete[] backupUpperBoundLevels;
        delete[] backupLowerBoundLevels;
        delete[] backupUpperBounds;
        delete[] backupLowerBounds;

        timeval end = Time::sampleMicro();
        _totalRestorationTimeMilli += Time::timePassed( start, end );

        printf( "\n\n\t\t !!! Restore tableau from backup DONE !!!\n" );
    }

    void adjustBasicVariables( const Set<unsigned> &shouldBeBasic, Set<unsigned> shouldntBeBasic, bool adjustAssignment = true )
    {
        unsigned count = 0;
        for ( const auto &entering : shouldBeBasic )
        {
            const Tableau::Entry *columnEntry = _tableau.getColumn( entering );
            const Tableau::Entry *current;

            bool done = false;
            while ( !done && columnEntry != NULL )
            {
                current = columnEntry;
                columnEntry = columnEntry->nextInColumn();

                unsigned leaving = current->getRow();
                if ( shouldntBeBasic.exists( leaving ) )
                {
                    double weight = FloatUtils::abs( getCell( leaving, entering ) );
                    if ( FloatUtils::lt( weight, NUMBERICAL_INSTABILITY_CONSTANT ) )
                    {
                        log( Stringf( "adjustBasicVariables: skipping a bad pivot: %.10lf\n",
                                      getCell( leaving, entering ) ) );
                        continue;
                    }

                    ++count;

                    done = true;
                    shouldntBeBasic.erase( leaving );

                    pivot( entering, leaving );
                    computeVariableStatus( leaving );

                    if ( adjustAssignment )
                    {
                        // leaving is now non-basic, so the invariants needs to be enforced.
                        if ( tooLow( leaving ) )
                            update( leaving, _lowerBounds[leaving].getBound() - _assignment[leaving], true );
                        else if ( tooHigh( leaving ) )
                            update( leaving, _upperBounds[leaving].getBound() - _assignment[leaving], true );

                        // If there's a Relu problem, fix leaving
                        if ( _reluPairs.isRelu( leaving ) )
                        {
                            unsigned b = _reluPairs.isB( leaving ) ? leaving : _reluPairs.toPartner( leaving );
                            unsigned f = _reluPairs.toPartner( b );

                            if ( ( !_dissolvedReluVariables.exists( f ) ) && reluPairIsBroken( b, f ) )
                            {
                                if ( ( !_basicVariables.exists( b ) ) && ( !_basicVariables.exists( f ) ) )
                                {
                                    if ( FloatUtils::isPositive( b ) )
                                        update( f, _assignment[b] - _assignment[f], true );
                                    else
                                        update( f, -_assignment[f], true );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void fixAllBrokenRelus()
    {
        for ( auto &pair : _reluPairs.getPairs() )
        {
            unsigned b = pair.getB();
            unsigned f = pair.getF();

            if ( ( !_dissolvedReluVariables.exists( f ) ) && reluPairIsBroken( b, f ) )
            {
                if ( ( !_basicVariables.exists( b ) ) && ( !_basicVariables.exists( f ) ) )
                {
                    if ( FloatUtils::isPositive( b ) )  // 让f与b保持一致
                        update( f, _assignment[b] - _assignment[f], true );
                    else
                        update( f, getLeakyValue() * _assignment[b] - _assignment[f], true );
                }
            }
        }
    }
    void fixAllBrokenRelus_temp()
    {
        for ( auto &pair : _reluPairs.getPairs() )
        {
            unsigned b = pair.getB();
            unsigned f = pair.getF();

            if ( ( !_dissolvedReluVariables.exists( f ) ) && reluPairIsBroken( b, f ) )
            {
                if ( ( !_basicVariables.exists( b ) ) && ( !_basicVariables.exists( f ) ) )
                {
                    if ( FloatUtils::isPositive( b ) )  // 让f与b保持一致
                        update( f, _assignment[b] - _assignment[f], true );
                    else
                        update( f, -_assignment[f], true );
                }
            }
        }
    }

    UseSlackVariables useSlackVariablesForRelus() const
    {
        return _useSlackVariablesForRelus;
    }

    Set<unsigned> getActiveRowSlacks() const
    {
        return _activeSlackRowVars;
    }

    Set<unsigned> getActiveColSlacks() const
    {
        return _activeSlackColVars;
    }

    double getSlackLowerBound( unsigned variable ) const
    {
        DEBUG({
                  if ( !_slackToLowerBound.exists( variable ) )
                  {
                      if ( ( ( _useSlackVariablesForRelus == USE_ROW_SLACK_VARIABLES ) &&
                             !_activeSlackRowVars.exists( variable ) )
                           ||
                           ( ( _useSlackVariablesForRelus == USE_ROW_AND_COL_SLACK_VARIABLES ) &&
                             !_activeSlackColVars.exists( variable ) ) )
                      {
                          printf( "Error! requested a slack lower bound on a non-slack variable (%u)!\n",
                                  variable );

                          exit( 1 );
                      }
                  }
              });

        return _slackToLowerBound.at( variable ).getBound();
    }

    double getSlackUpperBound( unsigned variable ) const
    {
        DEBUG({
                  if ( !_slackToUpperBound.exists( variable ) )
                  {
                      if ( ( ( _useSlackVariablesForRelus == USE_ROW_SLACK_VARIABLES ) &&
                             !_activeSlackRowVars.exists( variable ) )
                           ||
                           ( ( _useSlackVariablesForRelus == USE_ROW_AND_COL_SLACK_VARIABLES ) &&
                             !_activeSlackColVars.exists( variable ) ) )
                      {
                          printf( "Error! requested a slack upper bound on a non-slack variable (%u)!\n",
                                  variable );
                          exit( 1 );
                      }
                  }
              });

        return _slackToUpperBound.at( variable ).getBound();
    }

    unsigned slackToB( unsigned slack ) const
    {
        return _slackRowVariableToB.at( slack );
    }

    unsigned slackToF( unsigned slack ) const
    {
        return _slackRowVariableToF.at( slack );
    }

    int fixRelusInGlpkAssignment( int n, int m, int nonBasicEncoding, const int *head, const char *flags )
    {
        ++_fixRelusInGlpkAssignmentInvoked;

        unsigned nonBasic = _currentGlpkWrapper->glpkEncodingToVariable( nonBasicEncoding );
        if ( !activeReluVariable( nonBasic ) )
            return 0;

        // Check if partner is non-basic.
        unsigned partner = _reluPairs.toPartner( nonBasic );
        unsigned partnerEncoding = _currentGlpkWrapper->variableToGlpkEncoding( partner );

        bool partnerIsNonBasic = false;
        int partnerIndex = 1;
        for ( ; partnerIndex <= n-m; ++partnerIndex )
        {
            if ( head[m + partnerIndex] == (int)partnerEncoding )
            {
                partnerIsNonBasic = true;
                break;
            }
        }

        if ( !partnerIsNonBasic )
            return 0;

        char currentBound = flags[nonBasicEncoding];
        char partnerBound = flags[partnerEncoding];

        if ( currentBound != partnerBound )
        {
            // The partner needs to have its bound flipped
            if ( !_reluUpdateFrequency.exists( partner ) )
                _reluUpdateFrequency[partner] = 0;
            ++_reluUpdateFrequency[partner];
            if ( _reluUpdateFrequency[partner] > 5 )
            {
                ++_fixRelusInGlpkAssignmentIgnore;
                return 0;
            }

            ++_fixRelusInGlpkAssignmentFixes;
            return partnerIndex;
        }

        return 0;
    }

    void conflictAnalysisCausedPop()
    {
        ++_conflictAnalysisCausedPop;
    }

    void quit()
    {
        _quit = true;
    }

    void addTimeEvalutingGlpkRows( unsigned time )
    {
        _totalTimeEvalutingGlpkRows += time;
    }

    //// add by lzs
    double getLeakyValue() const {
        return _leakyValue;
    }
    void setLeakyValue(double leakyValue){
        _leakyValue = leakyValue;
    }
    //// add end

private:
    unsigned _numVariables;
    String _reluplexName;
    char *_finalOutputFile;
    FinalStatus _finalStatus;
    bool _wasInitialized;
    Tableau _tableau;
    Tableau _preprocessedTableau;
    VariableBound *_upperBounds;	//数组，存储所有变量的相应值
    VariableBound *_lowerBounds;	//数组，存储所有变量的相应值
    VariableBound *_preprocessedUpperBounds;	//数组，存储所有变量的相应值
    VariableBound *_preprocessedLowerBounds;	//数组，存储所有变量的相应值
    double *_preprocessedAssignment;	//数组，存储所有变量的相应值
    Set<unsigned> _basicVariables;	// Set
    Set<unsigned> _preprocessedBasicVariables;
    Map<unsigned, String> _variableNames;
    ReluPairs _reluPairs;
    SmtCore _smtCore;
    bool _useApproximations;
    bool _findAllPivotCandidates;
    unsigned _conflictAnalysisCausedPop;

    GlpkWrapper::GlpkAnswer _previousGlpkAnswer;

    bool _logging;
    bool _dumpStates;

    unsigned _numCallsToProgress;	//初始化为0，表示调用process的次数
    unsigned _numPivots;
    unsigned long long _totalPivotTimeMilli;
    unsigned long long _totalDegradationCheckingTimeMilli;
    unsigned long long _totalRestorationTimeMilli;
    unsigned long long _totalPivotCalculationCount;
    unsigned long long _totalNumBrokenRelues;
    unsigned _brokenRelusFixed;
    unsigned _brokenReluFixByUpdate;
    unsigned _brokenReluFixByPivot;
    unsigned _brokenReluFixB;
    unsigned _brokenReluFixF;
    unsigned _numEliminatedVars;
    unsigned _varsWithInfiniteBounds;
    unsigned _numStackSplits;
    unsigned _numStackMerges;
    unsigned _numStackPops;
    unsigned _numStackVisitedStates;
    unsigned _currentStackDepth;
    unsigned _minStackSecondPhase;
    unsigned _maximalStackDepth;
    unsigned long long _boundsTightendByTightenAllBounds;

    unsigned _almostBrokenReluPairCount;
    unsigned _almostBrokenReluPairFixedCount;

    unsigned _numBoundsDerivedThroughGlpk;
    unsigned _numBoundsDerivedThroughGlpkOnSlacks;
    unsigned long long _totalTightenAllBoundsTime;

    bool _eliminateAlmostBrokenRelus;

    Map<unsigned, VariableStatus> _varToStatus;

    Map<unsigned, ReluDissolutionType> _dissolvedReluVariables;     //
    Map<unsigned, ReluDissolutionType> _preprocessedDissolvedRelus; //

    bool _printAssignment;
    Set<unsigned> _eliminatedVars;  // 存储被消除的辅助变量，

    unsigned _numOutOfBoundFixes;
    unsigned _numOutOfBoundFixesViaBland;

    bool _useDegradationChecking;
    unsigned _numLpSolverInvocations;
    unsigned _numLpSolverFoundSolution;
    unsigned _numLpSolverNoSolution;
    unsigned _numLpSolverFailed;
    unsigned _numLpSolverIncorrectAssignment;
    unsigned long long _totalLpSolverTimeMilli;
    unsigned long long _totalLpExtractionTime;
    unsigned _totalLpPivots;
    unsigned _maxLpSolverTimeMilli;

    unsigned _numberOfRestorations;
    double _maxDegradation;

    unsigned long long _totalProgressTimeMilli;
    unsigned long long _timeTighteningGlpkBoundsMilli;

    GlpkWrapper *_currentGlpkWrapper;

    unsigned _relusDissolvedByGlpkBounds;

    Map<unsigned, VariableBound> _glpkStoredUpperBounds;
    Map<unsigned, VariableBound> _glpkStoredLowerBounds;

    double _glpkSoi;

    unsigned long long _storeGlpkBoundTighteningCalls;
    unsigned long long _storeGlpkBoundTighteningCallsOnSlacks;
    unsigned long long _storeGlpkBoundTighteningIgnored;

    unsigned _maxBrokenReluAfterGlpk;
    unsigned _totalBrokenReluAfterGlpk;
    unsigned _totalBrokenNonBasicReluAfterGlpk;

    UseSlackVariables _useSlackVariablesForRelus;
    Set<unsigned> _activeSlackRowVars;
    Set<unsigned> _activeSlackColVars;

    Map<unsigned, unsigned> _fToSlackRowVar;
    Map<unsigned, unsigned> _fToSlackColVar;

    Map<unsigned, unsigned> _slackRowVariableToF;
    Map<unsigned, unsigned> _slackRowVariableToB;
    Map<unsigned, VariableBound> _slackToLowerBound;
    Map<unsigned, VariableBound> _slackToUpperBound;

    Map<unsigned, unsigned> _reluUpdateFrequency;

    unsigned long long _fixRelusInGlpkAssignmentFixes;
    unsigned long long _fixRelusInGlpkAssignmentInvoked;
    unsigned long long _fixRelusInGlpkAssignmentIgnore;

    bool _maximalGlpkBoundTightening;
    bool _useConflictAnalysis;
    bool _temporarilyDontUseSlacks;

    bool _quit;
    bool _fullTightenAllBounds;
    bool _glpkExtractJustBasics;

    unsigned long long _totalTimeEvalutingGlpkRows;
    unsigned _consecutiveGlpkFailureCount;

    bool alreadySAT;     // add by lzs
    double _leakyValue;  // add by lzs


public:
    //// add by lzs
    void printTableauRow(unsigned row){
        _tableau.printRow(row);
    }
    //// add end

    void checkInvariants() const
    {
#ifndef DEBUG_ON
        return;
#endif

        // Table is in tableau form

        /*********check Tableau中是否有非法值，例如：
            1、basic值在Tableau中_column数组中对应位置一定不为null
            2、初始化Cell时最后一行是否满足initializeCell( 6, 6, -1.0 )的形式
            3、basic是否出现在其他basic的等式右边
        ***************/

        // _basicVariables 是 Set<unsigned>类型，&basic是每一个basic变量的地址
        for ( const auto &basic : _basicVariables )
        {
            // 判断Tableau中basic对应的位置是否有指向，如果有，则表示激活状态，如果为Null，则表示没有激活
            // 如果没有激活，则退出
            if ( !_tableau.activeColumn( basic ) )
            {
                printf( "Error: basic variable's column should be active! (var: %s)\n",
                        toName( basic ).ascii() );
                exit( 1 );
            }

            // 设置在Tableau的方阵中从basic指向的那一列为入口，即双向链表的Head,
            // importent:这里就要求在intialCell中构造的时候，最后一行必须是_reluplex->initializeCell( 6, 6, -1.0 )、_reluplex->initializeCell( 7, 7, -1.0 );的形式
            const Tableau::Entry *columnEntry = _tableau.getColumn( basic );

            // 对于basic变量，它只能出现在等式的左边，不能出现在右边，（因为但凡出现在右边的都会被代换掉
            // 所以，basic的列size必须为1，且head Entry对应的value必须是 -1.0 (构造时的设置必须与这里一致)
            // 另又由于构造要求，所以这个head Entry的row 和 column都是同一个值，所以要求columnEntry->getRow() == basic
            // 如果不等于，那么也要报错退出
            if ( ( _tableau.getColumnSize( basic ) != 1 ) ||
                 ( columnEntry->getRow() != basic ) ||
                 ( FloatUtils::areDisequal( columnEntry->getValue(), -1.0 ) ) )
            {
                printf( "Error: basic variable's column isn't right! (var: %s)\n",
                        toName( basic ).ascii() );
                printf( "Column size = %u\n", _tableau.getColumnSize( basic ) );

                exit( 1 );
            }

            const Tableau::Entry *rowEntry = _tableau.getRow( basic );
            const Tableau::Entry *current;

            while ( rowEntry != NULL )
            {
                current = rowEntry;
                rowEntry = rowEntry->nextInRow();

                if ( current->getColumn() == basic )
                    continue;

                // 判断当前basic的Tableau的row双向链表中，遍历每一行存在的Entry，取得他们的column，也就是等式右边存在的variable,
                // 如果发现也存在于basic中，报错退出
                if ( _basicVariables.exists( current->getColumn() ) )
                {
                    printf( "Error: a basic variable appears in another basic variable's row\n" );
                    exit( 1 );
                }
            }
        }

        /***** 检查 non-basic 要在上下界之内*****/
        // All non-basic are within bounds
        for ( unsigned i = 0; i < _numVariables; ++i )
        {
            if ( _varToStatus.get( i ) == ABOVE_UB || _varToStatus.get( i ) == BELOW_LB )
            {
                // Only basic variables can be out-of-bounds
                // 如果是non-basic，但有越界行为，则打印错误并退出
                if ( !_basicVariables.exists( i ) )
                {
                    printf( "Error: variable is out-of-bounds but is not basic! "
                            "(var: %s, value = %.10lf, range = [%.10lf, %.10lf])\n",
                            toName( i ).ascii(),
                            _assignment[i],
                            _lowerBounds[i].getBound(),
                            _upperBounds[i].getBound() );
                    if ( _reluPairs.isRelu( i ) )
                    {
                        unsigned partner = _reluPairs.toPartner( i );

                        printf( "This is also a relu variable (%s)\n", _reluPairs.isF( i ) ? "F" : "B" );
                        printf( "Relu partner is: %s. value = %.10lf, range = [%.10lf, %.10lf])\n",
                                toName( partner ).ascii(),
                                _assignment[partner],
                                _lowerBounds[partner].getBound(),
                                _upperBounds[partner].getBound() );
                    }

                    exit( 1 );
                }
            }
        }

        //
        // All relus that were split or merged have correct bounds and tableau shape
        for ( const auto &dissolved : _dissolvedReluVariables )
        {
            unsigned f = dissolved.first;
            unsigned b = _reluPairs.toPartner( f );

            double bUpper = _upperBounds[b].getBound();

            double fLower = _lowerBounds[f].getBound();
            double fUpper = _upperBounds[f].getBound();

            if ( dissolved.second == TYPE_SPLIT )
            {
                // F should be fixed, B's upper bound should be non-positive
                if ( !FloatUtils::isZero( fUpper ) || !FloatUtils::isZero( fLower ) )
                {
                    printf( "Error! After a split, F is not fixed at zero.\n"
                            "f = %u. Lower = %.15lf. Upper = %.15lf\n",
                            f, fLower, fUpper );
                    exit( 1 );
                }

                if ( FloatUtils::isPositive( bUpper ) )
                {
                    printf( "Error! After a split, B's upper bound is positive.\n"
                            "b = %u. Upper = %.15lf\n",
                            b, bUpper );
                    exit( 1 );
                }
            }
            else
            {
                // B should be eliminated, F should (naturally) have a non-negative lower bound
                if ( _tableau.getColumnSize( b ) != 0 )
                {
                    printf( "Error! After a merge, b's column is not of size 0! b = %u\n", b );
                    exit( 1 );
                }

                if ( FloatUtils::isNegative( fLower ) )
                {
                    printf( "Error! After a merge, F's lower bound is negative.\n"
                            "f = %u. Lower = %.15lf\n",
                            f, fLower );
                    exit( 1 );
                }
            }
        }
    }

    void printColumn( unsigned index )
    {
        printf( "\n\nDumping column for %s:\n", toName( index ).ascii() );

        const Tableau::Entry *columnEntry = _tableau.getColumn( index );
        while ( columnEntry != NULL )
        {
            printf( "\t<%u, %.5lf>\n", columnEntry->getRow(), columnEntry->getValue() );
            columnEntry = columnEntry->nextInColumn();
        }
    }

    static String statusToString( VariableStatus status )
    {
        switch ( status )
        {
            case ABOVE_UB: return "Above UB";
            case AT_UB: return "At UB";
            case BETWEEN: return "Between";
            case FIXED: return "Fixed";
            case AT_LB: return "At LB";
            case BELOW_LB: return "Below LB";
        }

        return "Unknown";
    }

    void tightenAllBounds()
    {
        log( "tightenAllBounds -- Starting\n" );

        timeval start = Time::sampleMicro();

        unsigned numLearnedBounds = 0;

        if ( !_fullTightenAllBounds )
        {
            Set<unsigned> copyOfBasics = _basicVariables;
            for ( const auto &basic : copyOfBasics )
            {
                if ( !_basicVariables.exists( basic ) )
                    continue;

                tightenBoundsOnRow( basic, numLearnedBounds );
            }
        }
        else
        {
            bool done = false;
            while ( !done )
            {
                bool needToRestart = false;
                Set<unsigned>::iterator basic = _basicVariables.begin();

                while ( ( basic != _basicVariables.end() ) && !needToRestart )
                {
                    needToRestart = tightenBoundsOnRow( *basic, numLearnedBounds );
                    ++basic;
                }

                if ( !needToRestart )
                    done = true;
            }
        }

        timeval end = Time::sampleMicro();
        _totalTightenAllBoundsTime += Time::timePassed( start, end );

        _boundsTightendByTightenAllBounds += numLearnedBounds;

        log( Stringf( "tightenAllBounds -- Done. Number of learned bounds: %u\n", numLearnedBounds ) );
    }

    bool tightenBoundsOnRow( unsigned basic, unsigned &numLearnedBounds )
    {
        const Tableau::Entry *row = _tableau.getRow( basic );
        const Tableau::Entry *tighteningVar;

        while ( row != NULL )
        {
            tighteningVar = row;

            row = row->nextInRow();

            double scale = -1.0 / tighteningVar->getValue();

            const Tableau::Entry *otherEntry = _tableau.getRow( basic );
            const Tableau::Entry *current;
            unsigned column;

            double max = 0.0;
            double min = 0.0;
            unsigned minBoundLevel = 0;
            unsigned maxBoundLevel = 0;

            while ( otherEntry != NULL )
            {
                current = otherEntry;
                otherEntry = otherEntry->nextInRow();

                column = current->getColumn();
                if ( column == tighteningVar->getColumn() )
                    continue;

                double coefficient = current->getValue() * scale;
                if ( FloatUtils::isPositive( coefficient ) )
                {
                    // Positive coefficient
                    min += _lowerBounds[column].getBound() * coefficient;
                    max += _upperBounds[column].getBound() * coefficient;

                    if ( _lowerBounds[column].getLevel() > minBoundLevel )
                        minBoundLevel = _lowerBounds[column].getLevel();
                    if ( _upperBounds[column].getLevel() > maxBoundLevel )
                        maxBoundLevel = _upperBounds[column].getLevel();
                }
                else
                {
                    // Negative coefficient
                    min += _upperBounds[column].getBound() * coefficient;
                    max += _lowerBounds[column].getBound() * coefficient;

                    if ( _lowerBounds[column].getLevel() > maxBoundLevel )
                        maxBoundLevel = _lowerBounds[column].getLevel();
                    if ( _upperBounds[column].getLevel() > minBoundLevel )
                        minBoundLevel = _upperBounds[column].getLevel();
                }
            }

            unsigned currentVar = tighteningVar->getColumn();

            if ( FloatUtils::lt( max, _upperBounds[currentVar].getBound() ) )
            {
                // Found an UB
                ++numLearnedBounds;
                updateUpperBound( currentVar, max, maxBoundLevel );
            }

            if ( FloatUtils::gt( min, _lowerBounds[currentVar].getBound() ) )
            {
                // Found a LB
                ++numLearnedBounds;
                // Tableau changed, need to restart
                if ( updateLowerBound( currentVar, min, minBoundLevel ) )
                    return true;
            }
        }

        // Don't need to restart
        return false;
    }

    void adjustGlpkAssignment( Map<unsigned, double> &assignment )
    {
        for ( auto &pair : assignment )
        {
            unsigned var = pair.first;
            double value = pair.second;

            if ( _basicVariables.exists( pair.first ) )
                continue;

            // Adjust variables to their bounds according to our precision
            if ( FloatUtils::gt( _lowerBounds[var].getBound(), value ) )
            {
                printf( "Adjust to lower bound. Var %u: value = %lf, bound = %lf\n",
                        var, value, _lowerBounds[var].getBound() );
                pair.second = _lowerBounds[var].getBound();
            }

            if ( FloatUtils::lt( _upperBounds[var].getBound(), value ) )
            {
                printf( "Adjust to upper bound. Var %u: value = %lf, bound = %lf\n",
                        var, value, _upperBounds[var].getBound() );
                pair.second = _upperBounds[var].getBound();
            }

            if ( ( pair.second != 0.0 ) && FloatUtils::isZero( pair.second ) )
            {
                pair.second = 0.0;
            }
        }
    }

    bool checkEquationsHold( unsigned basic, Map<unsigned, double> &assignment )
    {
        double result = 0.0;

        const Tableau::Entry *rowEntry = _tableau.getRow( basic );
        const Tableau::Entry *current;

        while ( rowEntry != NULL )
        {
            current = rowEntry;
            rowEntry = rowEntry->nextInRow();

            unsigned column = current->getColumn();

            if ( column != basic )
                result += assignment[column] * current->getValue();
            else
            if ( FloatUtils::areDisequal( current->getValue(), -1.0 ) )
            {
                printf( "Error! Basic's coefficient is not -1. It is: %lf\n", current->getValue() );
                exit( 2 );
            }
        }

        if ( !FloatUtils::areEqual( assignment[basic], result, GLPK_IMPRECISION_TOLERANCE ) )
        {
            printf( "Error! Mismatch between glpk answer and calculation for basic var = %s. "
                    "Calculated: %.10lf. Glpk: %.10lf.\n",
                    toName( basic ).ascii(), result, assignment[basic] );
            return false;
        }

        if ( FloatUtils::isZero( result ) )
            result = 0.0;

        assignment[basic] = result;
        return true;
    }

    void calculateBasicVariableValues()
    {
        for ( unsigned basic : _basicVariables )
            calculateBasicVariableValue( basic );
    }

    void calculateBasicVariableValue( unsigned basic )
    {
        double result = 0.0;

        const Tableau::Entry *rowEntry = _tableau.getRow( basic );
        const Tableau::Entry *current;

        while ( rowEntry != NULL )
        {
            current = rowEntry;
            rowEntry = rowEntry->nextInRow();

            unsigned column = current->getColumn();

            if ( column != basic )
                result += _assignment[column] * current->getValue();
        }

        if ( FloatUtils::isZero( result ) )
            result = 0.0;

        _assignment[basic] = result;
        computeVariableStatus( basic );
    }

    double checkDegradation()
    {
        timeval start = Time::sampleMicro();

        double max = 0.0;
        for ( unsigned basic : _preprocessedBasicVariables )
        {
            double degradation = checkDegradation( basic );
            if ( FloatUtils::gt( degradation, max ) )
                max = degradation;
        }

        if ( max > _maxDegradation )
            _maxDegradation = max;

        timeval end = Time::sampleMicro();
        _totalDegradationCheckingTimeMilli += Time::timePassed( start, end );

        return max;
    }

    double checkDegradation( unsigned variable )
    {
        double result = 0.0;

        const Tableau::Entry *rowEntry = _preprocessedTableau.getRow( variable );
        const Tableau::Entry *current;

        while ( rowEntry != NULL )
        {
            current = rowEntry;
            rowEntry = rowEntry->nextInRow();

            unsigned column = current->getColumn();

            if ( column != variable )
            {
                unsigned adjustedColumn = column;

                // If column is the b variable of a merged pair, take f instead
                if ( _reluPairs.isRelu( column ) && _reluPairs.isB( column ) )
                {
                    if ( _tableau.getColumnSize( column ) == 0 )
                        adjustedColumn = _reluPairs.toPartner( column );
                }

                result += _assignment[adjustedColumn] * current->getValue();
            }
        }

        unsigned adjustedVariable = variable;
        if ( _reluPairs.isRelu( variable ) && _reluPairs.isB( variable ) )
        {
            if ( _tableau.getColumnSize( variable ) == 0 )
                adjustedVariable = _reluPairs.toPartner( variable );
        }

        return FloatUtils::abs( result - _assignment[adjustedVariable] );
    }
};

void boundCalculationHook( int n, int m, int *head, int leavingBasic, int enteringNonBasic, double *basicRow )
{
    activeReluplex->storeGlpkBoundTightening( n, m, head, leavingBasic, enteringNonBasic, basicRow );
}

void iterationCountCallback( int count )
{
    activeReluplex->glpkIterationCountCallback( count );
}

void reportSoiCallback( double soi )
{
    activeReluplex->glpkReportSoi( soi );
}

int makeReluAdjustmentsCallback( int n, int m, int nonBasicEncoding, const int *head, const char *flags )
{
    return activeReluplex->fixRelusInGlpkAssignment( n, m, nonBasicEncoding, head, flags );
}

#endif // __Reluplex_h__

//
// Local Variables:
// compile-command: "make -C . "
// tags-file-name: "./TAGS"
// c-basic-offset: 4
// End:
//
