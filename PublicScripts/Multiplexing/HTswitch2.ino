
//int timer = 1000;
int seq[] = {1,1,1,0,1,0,0};
int seqLen = 7;
int cycleTime = 5000;
int injTime = 1000;
int waitTime = cycleTime - injTime;
//const int switchOut = A1;
//const int startOut = A0;

void setup() {
  // initialize digital pin LED_BUILTIN as an output.
  pinMode(LED_BUILTIN, OUTPUT);
  pinMode(2, OUTPUT);
  pinMode(3, OUTPUT);
}

void loop() {
  // put your main code here, to run repeatedly:
  unsigned long currentMillis = millis();
  unsigned long previousMillis = currentMillis;
  for (int i; i < seqLen; i++){
    if (i == 0){
      digitalWrite(3, HIGH);
      //if (currentMillis - previousMillis >= 1000) {
      //    digitalWrite(3, LOW);
      //}
    }

      /*
      analogWrite(startOut, 3);
      delay(100);
      analogWrite(startOut, 0);*/
      
    if (seq[i]==1){
      digitalWrite(LED_BUILTIN, HIGH);
      digitalWrite(2, HIGH); 
      delay(injTime);
      digitalWrite(LED_BUILTIN, LOW); 
      digitalWrite(2, LOW); 
      delay(waitTime);
    }
    else {
      delay(cycleTime);
    }
    if (i == 0){
      digitalWrite(3, LOW);
    }

  }  
}
