// If the collator should expect RDMA, then it has to enable the right switches
#ifdef COLLATOR_ENABLE_RDMA
#define ENABLE_RDMA_HANDSHAKES
#include "cmidirect.h"
#endif
#include "RDMAMessages.h"
#include "charm++.h"
#include <string>
#include <sstream>

#ifndef MESSAGE_DATA_COLLATOR_H
#define MESSAGE_DATA_COLLATOR_H

/*
 * @addtogroup PairCalculator
 *    @{
 */

namespace cp {
  namespace paircalc {

#ifdef DEBUG_MESSAGEDATACOLLATOR_ALL
#define DEBUG_MESSAGEDATACOLLATOR
#define DEBUG_MESSAGEDATACOLLATOR_CREATION
#define DEBUG_MESSAGEDATACOLLATOR_RDMA
#define DEBUG_MESSAGEDATACOLLATOR_EACHARRIVAL
#endif

    /** Class that buffers incoming data (via messages/RDMA) till it counts a pre-specified number of arrivals, and then spits the collated data out 
     * (packaged in the same incoming message type) via a CkCallback. 
     *
     * Assumptions:
     *    - All messages carry the same amount of data
     *    - Each batch of data arrives completely only via RDMA or msgs (no interleaving of different methods within a single batch of messages)
     *    - Incoming messages are not deleted after processing. This class is not a chare and will let the user code deal with such formalities as needed
     *    - If this class manages the data buffer, then the post-collation callback goes to a chare on the same PE so that the output msg is not deleted/moved/resized in any way
     *    - The first batch of messages immediately after an RDMA setup (which will involve a bunch of setupRDMA() calls) WILL be via RDMA (and not messages)
     *    - expectNext() will be called before every new batch of RDMA arrivals (and only after the current batch is done)
     *    - Every time the amount of data in an imminent batch will be more than in any previous batch, the buffer is reallocated. A fresh RDMA setup is required for subsequent RDMA arrivals
     *    - Incoming data is treated as 2d (nRows x nCols). This is not binding. Higher / Lower dimensional data can be sent in as a contiguous array as long as msgType has the same API
     *    - Some effort has been made to tolerate different data sizes across different batches, but no testing / guarantees as yet
     *
     *
     * @note: Class msgType should provide the following methods:
     * - int numRows() const       : returns the number of rows in the dataType array in the message
     * - int numCols() const       : returns the number of columns in the dataType array in the message
     * - int sender() const        : returns the sender ID or a representation of the senderID as an int
     * - dataType* data()          : returns a pointer to the data held by the message. Does not strictly
     *                               have to be of dataType as it is only used by CmiMemcpy. 
     * - msgType(const int senderID, const int nCols, const int nRows) : Used to create the output msg that holds the collated data
     */
    template <class msgType, typename dataType>
                                               class MessageDataCollator 
    {
      public:
        // ----------------------------------------- Cons/Des-truction -----------------------------------------
        /// @warning: Do not use this. This exists only to appease charm migration code which screams for default constructors.
        MessageDataCollator() {}
        /** Collators should be constructed thus. Copying has not been restricted so that collator classes can be passed in as
         * arguments. Hence, be careful when copying collators to make sure multiple collators dont use the same memory region as
         * a data buffer and kill you.
         */
        MessageDataCollator(std::string debugID, CkCallback callbk, const unsigned int nExpected, bool manageBuffer=false, const unsigned int offsetStart=0):
          myName(debugID),
          uponCompletion(callbk),
          isBufferMine(manageBuffer),
          outputMsg(0),
          dataBuffer(0),
          bufferSizeInUnits(0),
          numUnitsInMsg(0),
          numExpected(nExpected),
          numReceived(0),
          numBatchesProcessed(0),
          isNextBatchViaRDMA(false),
          isThisBatchViaRDMA(false),
          minSender(offsetStart),
          sampleMsg(0)    
      {
#ifdef DEBUG_MESSAGEDATACOLLATOR_CREATION
        CkPrintf("MessageDataCollator%s instantiated: numExpected = %d and buffer management = %d on(1)/off(0)\n",myName.c_str(),numExpected,isBufferMine);
#endif
#ifdef COLLATOR_ENABLE_RDMA
        /** If RDMA is enabled, the client code CANNOT own the buffer and do what it pleases with it. The buffer will 
         * have to be managed by me. Hence, this check. 
         * As an aside, for RDMA, I will choose to never delete/ reallocate the buffer as that will break RDMA. A repurcussion is that
         * setupRDMA() has to be called everytime the message size changes as that will prompt me to delete and reallocate the buffer
         */    
        if (!isBufferMine)
          CkAbort("MessageDataCollator: You cannot enable RDMA and want to manage the data buffer too. The buffer will have to be mine.");
        /// Reserve enough space for the RDMA handles that will be created for each data sender 
        RDMAhandles.reserve(numExpected); 
#endif
      }

        /// Destructor
        ~MessageDataCollator()
        {
          if (isBufferMine)
            delete outputMsg;
          delete sampleMsg;
        }

        /// PUP method @todo: Define a PUP method
        // void pup(PUP::er &p); 

        // ---------------------------------------- Message Handling ----------------------------------------
        /// Makes this a functor so that it can be fed messages by any entity that just knows the MessageHandler API. 
        void operator() (msgType *msg);

        // ---------------------------------------- RDMA support ----------------------------------------
        /// Called by a data sender to setup an RDMA link with me (data receiver)
        template<class tokenType>
          void setupRDMA(RDMASetupRequestMsg<tokenType> *msg, tokenType *replyToken=0);
        /// RDMA equivalent of msg handling. When data arrives via RDMA, this is called for handling it. Gets called from the CmiDirect code via a callback. 
        void operator() (void);
        /// Fed to CmiDirect as a callback to be triggered when some data is deposited by RDMA
        static void collatorCallback(void *thisObj);
        /// Ask me to notify converse that we're ready for the next batch of messages.
        void expectNext();

        // ---------------------------------------- Utility get methods ----------------------------------------
        /// Get the number of messages this collator has been configured to buffer
        inline unsigned int getNumExpected() const { return numExpected; }
        /// Check the number of messages that have arrived
        inline unsigned int getNumReceived() const { return numReceived; }
        /// Get the number of batches of incoming messages already processed
        inline unsigned int getNumBatchesDone() const { return numBatchesProcessed; }
        /// Get a copy of a sample message from the last completed batch of messages. Copy managed by the user
        inline msgType* getSampleMsg(); 

      private:
        /// A string that encodes some kind of ID supplied by the owner of an instance
        std::string myName;
        /// Callback that is passed in to trigger the post-collation work. 
        CkCallback uponCompletion;
        /// If buffer is mine, I manage it, viz. reuse the same buffer for every batch of msgs where possible. If not, I allocate a fresh buffer and let the client manage it.
        bool isBufferMine;
        /// The message (which has the data buffer) in which the incoming data is collated
        msgType *outputMsg;
        /// A pointer to the memory region (contained in outputMsg) where incoming data is collated
        dataType *dataBuffer;
        /// Handle to a sample message copied from every batch of incoming messages
        msgType *sampleMsg;
        /// Counters to keep track of the messages
        unsigned int numExpected, numReceived, numBatchesProcessed;
        /// Variables related to the message/buffer sizes 
        unsigned int numRows, numCols, numUnitsInMsg, bufferSizeInUnits;
        /** Senders should have an ID accessed by msgType::sender(). This stores the minimum value of sender(), 
         * so that offsets into the memory buffer can be computed
         *  
         * @todo: A better way to accomplish that leaves the collator class less rigid might be better
         */
        unsigned int minSender;
        /// A vector of RDMA handles, one for each sender
        CkVec<rdmaHandleType> RDMAhandles;
        /// Booleans indicating if the next/current batch of messages are via RDMA
        bool isNextBatchViaRDMA, isThisBatchViaRDMA;
    };



    template <class msgType,typename dataType>
      void MessageDataCollator<msgType,dataType>::operator() (msgType *msg)
      {
        /** @todo: Debate if some kind of check is needed to ensure that the message received is the correct iteration etc.
         * For e.g, if the message defined a uniqID() method that would return a unique tag indicating the group the message
         * belonged to, then all messages that belonged to this group could be buffered together. uniqID() could return
         * iteration number if messages from a single iteration are to be grouped, or the sender ID if messages from a single
         * sender are to be grouped together. Other possible scenarios where it is possible to receive messages from separate
         * collation groups might occur.
         */

        /// If we're in the middle of a batch of messages
        if (numReceived != 0)
          CkAssert(isThisBatchViaRDMA == false);
        /// else, if this is a fresh batch of messages
        else 
        {
          /// First, ensure that I was not expecting the next batch of data via RDMA 
          CkAssert(isNextBatchViaRDMA == false);
          /// Since, this is the first message in a new batch, flag this whole batch as NOT via RDMA
          isThisBatchViaRDMA = false;
          /// Get the message size
          numRows       = msg->numRows();
          numCols       = msg->numCols();
          numUnitsInMsg = numRows * numCols; 
          /// Get the required buffer size
          int reqdBufferSizeInUnits = numExpected * numUnitsInMsg;
          /// If buffer exists AND I manage it AND its not big enough for this batch of messages, delete it
          if (dataBuffer && isBufferMine && reqdBufferSizeInUnits > bufferSizeInUnits)
          {
#ifdef DEBUG_MESSAGEDATACOLLATOR
            CkPrintf("MessageDataCollator%s: Deleting buffer as I need space for %d units. Available = %d units\n",myName.c_str(),reqdBufferSizeInUnits,bufferSizeInUnits);
#endif
            delete outputMsg;
            outputMsg = 0; 
            dataBuffer = 0;
          }
          /// If the buffer is unallocated ... 
          if (!dataBuffer)
          {
            /// Get some memory to hold the message data
            bufferSizeInUnits =  reqdBufferSizeInUnits;
#ifdef DEBUG_MESSAGEDATACOLLATOR
            CkPrintf("MessageDataCollator%s: Allocating buffer of size %d units for %d data packets each having %d units in %d rows x %d cols.\n",myName.c_str(),bufferSizeInUnits,numExpected,numUnitsInMsg,numRows,numCols);
#endif
            outputMsg  = new (numExpected*numRows*numCols) msgType(-1,numCols,numRows*numExpected);
            dataBuffer = outputMsg->data();
            /// Zero out the memory region
            bzero(dataBuffer,bufferSizeInUnits * sizeof(dataType));
          } /// endif
        }

        /// Ensure that the output message and its data buffer exist
        CkAssert( (outputMsg != NULL) && (dataBuffer != NULL) );
        /// Ensure that this message has the same size as the others in this batch
        CkAssert( (numRows == msg->numRows()) && (numCols == msg->numCols()) );

        /// Compute the offset for this message. This will depend on the sender index and thisIndex. 
        int offset = msg->sender() - minSender;
        /// Copy the data from the message into the contiguous data buffer
        CmiMemcpy(&(dataBuffer[offset*numUnitsInMsg]), msg->data(), numUnitsInMsg *sizeof(dataType));

        /// Increment the tally of received messages 
        numReceived++;

        /// If the expected number have arrived ...
        if (numReceived >= numExpected)
        {
#ifdef DEBUG_MESSAGEDATACOLLATOR
          CkPrintf("MessageDataCollator%s: Collated data from %d of %d messages at %p. Gonna trigger my client\n",myName.c_str(),numReceived,numExpected,dataBuffer);
#endif
          /// Increment the number of message batches that have been processed
          numBatchesProcessed++;
          /// Create a copy of this message for a sample from this batch
          if (sampleMsg)
            delete sampleMsg;
          sampleMsg = reinterpret_cast<msgType*> ( CkCopyMsg(reinterpret_cast<void**>(&msg)) );
          /// Call the trigger functor
          uponCompletion.send(outputMsg);
          /// Reset the number received so that we can start collating all over again
          numReceived = 0;
          /// If you are not managing the buffer, then you need to forget about it
          if (!isBufferMine)
          {
            outputMsg  = 0;
            dataBuffer = 0;
          }
        } /// endif
        else
        {
#ifdef DEBUG_MESSAGEDATACOLLATOR_EACHARRIVAL
          CkPrintf("MessageDataCollator%s: Received message %d of %d.\n",myName.c_str(),numReceived,numExpected);
#endif
        }
        /// Do NOT delete the msg; MessageDataCollator is not a chare and doesnt worry about the message formalities. Let the user code handle this
      }



    /** Useful in case the msg has info other than the data() field. These can then 
     * be extracted in a case-specific manner, by the client code directly. 
     * 
     * @note: This class makes no guarantees on which message will be held as a sample msg.
     */
    template <class msgType,typename dataType>
      inline msgType* MessageDataCollator<msgType,dataType>::getSampleMsg() 
      { 
        if (sampleMsg) 
          return reinterpret_cast<msgType*> ( CkCopyMsg(reinterpret_cast<void**>(&sampleMsg)) );
        else
          return 0; 
      }



    template <class msgType,typename dataType> template <class tokenType>
      void MessageDataCollator<msgType,dataType>::setupRDMA(RDMASetupRequestMsg<tokenType> *msg, tokenType *replyToken)
      {
#ifndef COLLATOR_ENABLE_RDMA
        CkAbort("MessageDataCollator aborting because someone called an RDMA method when it has been turned off for me...\n");
#else
#ifdef DEBUG_MESSAGEDATACOLLATOR_RDMA
        std::stringstream dbgStr; 
        dbgStr<<"MessageDataCollator"<<myName<<": Received RDMA setup request from sender "<<msg->sender()<<" on proc "<<msg->senderPE()<<"."
          <<"\n\tRDMA token at Receiver: "<<*replyToken<<" dataSize = "<<msg->numDataUnits();
        CkPrintf("%s\n",dbgStr.str().c_str());
#endif
        /// You usually will not want to setup an RDMA link in the middle of receiving a batch of messages
        CkAssert(numReceived == 0);
        /// Flag the immediate next batch of messages as via RDMA ( which I expect, the user will respect :) ) 
        isNextBatchViaRDMA = true; 

        /// Get the message size
        numCols = msg->numDataUnits();
        numRows = 1;
        numUnitsInMsg = numRows * numCols;
        /// Get the required buffer size
        int reqdBufferSizeInUnits = numExpected * numUnitsInMsg;
        /// If buffer exists AND I manage it AND its not big enough for this batch of data, delete it
        if (dataBuffer && isBufferMine && reqdBufferSizeInUnits > bufferSizeInUnits)
        {
#ifdef DEBUG_MESSAGEDATACOLLATOR_RDMA
          CkPrintf("MessageDataCollator%s: Deleting buffer as I need space for %d units. Available = %d units\n",myName.c_str(),reqdBufferSizeInUnits,bufferSizeInUnits);
#endif
          delete outputMsg;
          outputMsg  = 0;
          dataBuffer = 0;
        }
        /// If the buffer is unallocated ... 
        if (!dataBuffer)
        {
          /// Get some memory to hold the data
          bufferSizeInUnits =  reqdBufferSizeInUnits;
#ifdef DEBUG_MESSAGEDATACOLLATOR_RDMA
          CkPrintf("MessageDataCollator%s: Allocating buffer of size %d units for %d data packets each having %d units in %d rows x %d cols.\n",myName.c_str(),bufferSizeInUnits,numExpected,numUnitsInMsg,numRows,numCols);
#endif
          outputMsg  = new (numExpected*numRows*numCols) msgType(-1,numCols,numRows*numExpected);
          dataBuffer = outputMsg->data();
          /// Zero out the memory region
          bzero(dataBuffer,bufferSizeInUnits * sizeof(dataType));
        } /// endif

        /// Compute the offset for this sender.This may depend on (for eg) the sender index and thisIndex
        int offset = msg->sender() - minSender;
        /// Get the region of the data buffer that will store data from this sender 
        dataType *recvLoc = &(dataBuffer[offset*numUnitsInMsg]);
        /// Create a function pointer to give converse for a callback
        void (*fnPtr)(void *) = MessageDataCollator<msgType,dataType>::collatorCallback;
        //double const actualQNaN = std::numeric_limits<double>::quiet_NaN();
        double const actualQNaN = 999999999999999999999999999999999999999999999999.9999999999;
        /// Create the RDMA handle for this sender-receiver pair
        rdmaHandleType rHandle = CmiDirect_createHandle(msg->senderPE(),recvLoc,numUnitsInMsg*sizeof(dataType),
            fnPtr,reinterpret_cast<void*>(this),actualQNaN);
        /// Store this handle for reference
        RDMAhandles.push_back(rHandle);
        /// Determine if the handshake token to be sent in the response has been supplied explicitly or should simply be copied from the request msg 
        tokenType rToken = msg->token();
        if (replyToken !=0)
          rToken = *replyToken;

        /// Prepare an RDMA request acceptance notification (and include the RDMA handle in the message for the data sender's use)
        RDMASetupConfirmationMsg<tokenType> *acceptanceMsg = new RDMASetupConfirmationMsg<tokenType>(rToken,rHandle);
        /// Notify the data sender of acceptance, to complete the RDMA setup handshake
        msg->senderCallback().send(acceptanceMsg);
        /// Do NOT delete the msg; MessageDataCollator is not a chare and doesnt worry about the message formalities. Let the user code handle this
#endif // ENABLE_COLLATOR_RDMA
      }



    template <class msgType,typename dataType>
      inline void MessageDataCollator<msgType,dataType>::operator() (void)
      {
#ifndef COLLATOR_ENABLE_RDMA
        CkAbort("MessageDataCollator aborting because someone called an RDMA method when it has been turned off for me...\n");
#else
        /// If this is a fresh batch of messages, flag it as via RDMA, else check if it is
        if (numReceived != 0)
          CkAssert(isThisBatchViaRDMA == true);
        else
        {
          /// Ensure that I was expecting this batch to be via RDMA (@todo: does this hold for the very first batch too?)
          CkAssert(isNextBatchViaRDMA == true);
          /// Flag this batch as via RDMA
          isThisBatchViaRDMA = true;
          /// Next batch is via messages unless I am explicitly told otherwise
          isNextBatchViaRDMA = false;
        }

        /// Increment the tally of completed RDMAs 
        numReceived++;

        /// If the expected number have arrived ...
        if (numReceived >= numExpected)
        {
#ifdef DEBUG_MESSAGEDATACOLLATOR_RDMA
          CkPrintf("MessageDataCollator%s: Collated %d of %d RDMA data packets. Gonna trigger my client\n",myName.c_str(),numReceived,numExpected);
#endif
          /// Increment the number of message batches that have been processed
          numBatchesProcessed++;
          /// Batches of RDMA cannot have sample messages. Delete if older one exists.
          if (sampleMsg)
            delete sampleMsg;
          sampleMsg = 0;
          /// Trigger the post-collation work via the callback
          uponCompletion.send(outputMsg);
          /// Reset the number received so that we can start collating all over again
          numReceived = 0;
        } /// endif
        else
        {
#ifdef DEBUG_MESSAGEDATACOLLATOR_EACHARRIVAL
          CkPrintf("MessageDataCollator%s: Someone deposited data via RDMA. Collated %d of %d.\n",myName.c_str(),numReceived,numExpected);
#endif
        }
#endif
      }



    template <class msgType,typename dataType>
      inline void MessageDataCollator<msgType,dataType>::collatorCallback(void *thisObj)
      {
#ifndef COLLATOR_ENABLE_RDMA
        CkAbort("MessageDataCollator aborting because someone called an RDMA method when it has been turned off for me...\n");
#else
        MessageDataCollator<msgType,dataType> *thisCollator = reinterpret_cast<MessageDataCollator<msgType,dataType>*>(thisObj);
        (*thisCollator)();
#endif
      }



    template <class msgType,typename dataType>
      inline void MessageDataCollator<msgType,dataType>::expectNext() 
      {
#ifndef COLLATOR_ENABLE_RDMA
        CkAbort("MessageDataCollator aborting because someone called an RDMA method when it has been turned off for me...\n");
#else
        /// Ensure that we're not in the middle of receiving a batch of messages @todo: Should this be a non-assert if condition thats always checked?
        CkAssert( numReceived==0 || numReceived>=numExpected );
        /// Set the flag the will tell me to not allow data via messages for the next batch
        isNextBatchViaRDMA = true; 
#ifdef DEBUG_MESSAGEDATACOLLATOR_RDMA
        CkPrintf("MessageDataCollator%s: Expecting the next batch of data via RDMA. Notifying the runtime to watch %d receive buffers.\n",myName.c_str(),RDMAhandles.size());
#endif
        /// Notify the runtime to expect data arrival via RDMA 
        for (int i=0;i<RDMAhandles.size();i++)
          CmiDirect_ready(&RDMAhandles[i]);
#endif
      }	


  } // end namespace paircalc
} // end namespace cp
/*@}*/
#endif // MESSAGE_DATA_COLLATOR_H

