/*******************************************************************************
  Main Source File

  Company:
    Microchip Technology Inc.

  File Name:
    main.c

  Summary:
    This file contains the "main" function for a project.

  Description:
    This file contains the "main" function for a project.  The
    "main" function calls the "SYS_Initialize" function to initialize the state
    machines of all modules in the system
 *******************************************************************************/

// *****************************************************************************
// *****************************************************************************
// Section: Included Files
// *****************************************************************************
// *****************************************************************************

#include <cstdio>
#include <cstdint>
#include <cstdarg>
#include <cstddef>                     // Defines NULL
#include <cstdbool>                    // Defines true
#include <cstdlib>                     // Defines EXIT_FAILURE
#include <cmath>
#include "definitions.h"                // SYS function prototypes
#include "edge-impulse-sdk/porting/ei_classifier_porting.h"
#include "edge-impulse-sdk/classifier/ei_run_classifier.h"
#include "edge-impulse-sdk/dsp/numpy.hpp"
#include "model-parameters/model_metadata.h"
#include "app.h"

#if EI_CLASSIFIER_SLICE_SIZE != MAX_AUDIO_NUM_SAMPLES
#error "MAX_AUDIO_NUM_SAMPLES ("#MAX_AUDIO_NUM_SAMPLES#") in app_config.h must match EI_CLASSIFIER_SLICE_SIZE ("#EI_CLASSIFIER_SLICE_SIZE#")"
#endif

/**
 * Wrapper around malloc
 */
void *ei_malloc(size_t size) {
    return malloc(size);
}

/**
 * Wrapper around calloc
 */
void *ei_calloc(size_t nitems, size_t size) {
    return calloc(nitems, size);
}

/**
 * Wrapper around free
 */

void ei_free(void *ptr) {
    free(ptr);
}

EI_IMPULSE_ERROR ei_run_impulse_check_canceled() {
    return EI_IMPULSE_OK;
}

EI_IMPULSE_ERROR ei_sleep(int32_t us) {
    TC2_Timer32bitCounterSet(0);
    while (TC2_Timer32bitCounterGet() < us) {};
    
    return EI_IMPULSE_OK;
}


uint64_t ei_read_timer_ms() {
    return TC2_Timer32bitCounterGet() / 1000;
}

uint64_t ei_read_timer_us() {
    return TC2_Timer32bitCounterGet() / 1;
}

void ei_printf(const char *format, ...) {
    va_list myargs;
    va_start(myargs, format);
    vprintf(format, myargs);
    va_end(myargs);
}

void ei_printf_float(float f) {
    ei_printf("%f", f);
}

extern "C" int get_feature_data(size_t offset, size_t length, float *out_ptr) {
    DRV_I2S_DATA *ptr = micBuffer[appData.rxBufferIdx];

    if (offset + length > EI_CLASSIFIER_SLICE_SIZE) {
        return EI_IMPULSE_CANCELED;
    }
    
    for (size_t i=0; i < length; i++) { 
        out_ptr[i] = (float) ptr[i + offset].rightData / 65536; // Convert Q31 to Q15
        
    }
    
    return EI_IMPULSE_OK;
}

extern "C" int _open (const char *buf, int flags, int mode)
{
    return 0;
}

void LEDTicker_Callback(uintptr_t context) {
    static int ticker = 0;
    int tickrate = *((int *) context);
    
    if (tickrate == 0) {
        ticker = 0;
        LED1_Off();
    }
    else if (++ticker == tickrate) {
        LED1_Toggle();
        ticker = 0;
    }
}
// *****************************************************************************
// *****************************************************************************
// Section: Main Entry Point
// *****************************************************************************
// *****************************************************************************


int main ( void )
{
    int tickrate = 0;
    
    /* Initialize all modules */
    SYS_Initialize ( NULL );
    
    /* Register and start the LED ticker */
    DRV_HANDLE tmrHandle;
    tmrHandle = SYS_TIME_CallbackRegisterMS(LEDTicker_Callback, 
                    (uintptr_t) &tickrate, 100, SYS_TIME_PERIODIC);
    
    TC2_TimerStart();
    
    ei::signal_t signal;
    signal.total_length = EI_CLASSIFIER_SLICE_SIZE;
    signal.get_data = &get_feature_data;
    
    ei_impulse_result_t result = { 0 };
    EI_IMPULSE_ERROR ei_status = EI_IMPULSE_OK;
    
    run_classifier_init();
    
    int votes[EI_CLASSIFIER_SLICES_PER_MODEL_WINDOW] = {0};
    int runningvotes = 0;
    while ( true )
    {
        /* Maintain state machines of all polled MPLAB Harmony modules. */
        SYS_Tasks ( );

        if (!appData.bufferReady) {
            continue;
        }

        ei_status = run_classifier_continuous(&signal, &result, false);
        if (ei_status != EI_IMPULSE_OK) {
            printf("run_classifier returned: %d\r\n", ei_status);
            break;
        }

        float maxval = 0;
        int maxidx = -1;
        for (size_t ix = 0; ix < EI_CLASSIFIER_LABEL_COUNT; ix++) {
            if (result.classification[ix].value > maxval) {
                maxval = result.classification[ix].value;
                maxidx = ix;
            }
        }        

        printf("Predictions (DSP: %d ms., Classification: %d ms., Anomaly: %d ms.): \r\n",
        result.timing.dsp, result.timing.classification, result.timing.anomaly);                
        for (size_t ix = 0; ix < EI_CLASSIFIER_LABEL_COUNT; ix++) {
            printf("    %s: %f\r\n", result.classification[ix].label, result.classification[ix].value);
        }

        // Roll vote buffer
        runningvotes -= votes[0];
        for (int i=1; i <  EI_CLASSIFIER_SLICES_PER_MODEL_WINDOW; i++)
            votes[i-1] = votes[i];
        
        // Threshold confidence score
        if (maxval > 0.7 && !strcmp(result.classification[maxidx].label, "vacuum"))
            votes[EI_CLASSIFIER_SLICES_PER_MODEL_WINDOW-1] = 1;
        else
            votes[EI_CLASSIFIER_SLICES_PER_MODEL_WINDOW-1] = 0;
        
        runningvotes += votes[EI_CLASSIFIER_SLICES_PER_MODEL_WINDOW-1];
        
        //printf("Votes: %d\r\n", runningvotes);
        if (runningvotes > EI_CLASSIFIER_SLICES_PER_MODEL_WINDOW/2) {
            tickrate = 1;
        }
        else {
            tickrate = 0;
        }
    }

    /* Execution should not come here during normal operation */

    return ( EXIT_FAILURE );
}


/*******************************************************************************
 End of File
*/

