
//int timer = 1000;

// standard 5 bit, 8 zp
//int seq[] = {0,0,0,0,1,0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,0,0,0,0,0,0};

// standard 5 bit, 10 zp
int seq[] = {0,0,0,0,1,0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0};

// standard 5 bit, 15 zp
//int seq[] = {0,0,0,0,1,0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

// inverse 5 bit, 10 zp
//int seq[] = {1,1,1,1,0,1,1,0,1,0,0,1,1,0,0,0,0,0,1,1,1,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0};  //{1,1,1,0,1,0,0,0}; {1,1,1,0,1,0,0,0,0,0,0,0,0,0,0};   {0,0,0,0,1,0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,0,0,0,0,0};      //{1,1,1,0,1,0,0,0,0,0};  //{0,0,0,0,1,0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1}; {1,1,1,0,1,0,0,0,0,0,0,0}   //{1,1,1,0,1,0,0};    //{0,0,0,0,1,0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1}

// 1 bit
//int seq[] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 

// 3 bit
//int seq[] = {1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0};

// 6 bit
//int seq[] = {0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,0,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,0,1,0,1,0,1,1,1,1,1,1}

// 6 bit, 15 zero padds
//int seq[] = {0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,0,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,0,1,0,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

//1,1,1,1,0,1,1,0,1,0,0,1,1,0,0,0,0,0,1,1,1,0,0,1,0,0,0,1,0,1,0
//0,0,0,0,1,0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1
int seqLen = 78;
int cycleTimeSeconds = 60.0L;
unsigned long cycleTime = 60000;
unsigned long injTime = 1500;
unsigned long waitTime = cycleTime - injTime;
int switchPin = 2;
int startPin = 5;
//const int switchOut = A1;
//const int startOut = A0;

void setup() {
  // initialize digital pin LED_BUILTIN as an output.
  pinMode(LED_BUILTIN, OUTPUT);
  pinMode(switchPin, OUTPUT);
  pinMode(startPin, OUTPUT);
}

void loop() {
  // put your main code here, to run repeatedly:
  unsigned long currentMillis = millis();
  unsigned long previousMillis = currentMillis;
  for (int i; i < seqLen; i++){
    if (i == 0){
      digitalWrite(startPin, HIGH);
    }

    if (seq[i]==1){
      digitalWrite(LED_BUILTIN, HIGH);
      digitalWrite(switchPin, HIGH); 
      delay(injTime);
      digitalWrite(LED_BUILTIN, LOW); 
      digitalWrite(switchPin, LOW); 
      delay(waitTime);
    }
    else {
      delay(cycleTime);
    }
    if (i == 0){
      digitalWrite(startPin, LOW);
    }

  }  
}
