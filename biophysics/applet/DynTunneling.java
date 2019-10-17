import java.awt.*;
import java.applet.*;

public class DynTunneling extends Applet implements Runnable {
  
  int NumberofStates = 5;
  double Steplength = 0.05;
  double[][] Hamiltonian = new double[NumberofStates][NumberofStates];
  Complex[] Psi = new Complex[NumberofStates];
  double Frequency[] = new double[NumberofStates];
  double Amplitude[] = new double[NumberofStates];
  double Occupancy[] = new double[NumberofStates];
  int PulseHeight[] = new int[NumberofStates];
  int HamilHeight[] = new int[NumberofStates];
  boolean printHamiltonian = false;
  boolean printInitialPsi = false;
  int outputInterval = 200;
  Thread myThread = null;
  int delay = 0;
  int AppletWidth,AppletHeight,xMargin,yMargin;
  int IndexofStep = 0;
  int IncrementofX, IncrementofY;
  boolean threadSuspended = false;
  double time = 0.0;
  double BarrierHeight=4.0;
  double Period = 100.0;
  boolean beginIt = false;
  Graphics offGraphics,onGraphics;
  Image offImage;

  public void init () {

    AppletWidth = size().width;
    AppletHeight = size().height;
    xMargin = AppletWidth / 10;
    yMargin = AppletHeight / 10;
    IncrementofX = (AppletWidth - 2 * xMargin)/(NumberofStates-1);
    IncrementofY = (AppletHeight -2 * yMargin);

    getHamiltonian();
    resetPsi();
    if(!beginIt) {
      Converter();
      delay=1000;
    }
    getVibrations();

    Panel panel = new Panel ();
    setLayout (new BorderLayout());
    panel.setLayout (new FlowLayout());
    add("North",panel);
    panel.add(new Button("Reset"));
    panel.add(new Button("Start"));
    Choice c = new Choice();
    c.addItem("Frequencies");
    c.addItem("Highest");
    c.addItem("Higher");
    c.addItem("High");
    c.addItem("Middle");
    c.addItem("Low");
    panel.add(c);

    show();
    onGraphics=getGraphics();
    update(onGraphics);
    offImage=this.createImage(AppletWidth,AppletHeight);
    offGraphics=offImage.getGraphics();

  }

  public void start () {
    if(myThread==null) {
      myThread = new Thread(this);
      myThread.start();
    }    
  }

  public void stop () {
    myThread = null;
  }

  public void run () {
    Thread.currentThread().setPriority(Thread.MIN_PRIORITY);
    while (myThread!=null) {
      if(beginIt) {RungeKutta();}
      if(IndexofStep%outputInterval==0) { 
        offGraphics.setColor(Color.black);
        offGraphics.fillRect(1,1,AppletWidth-2,AppletHeight-2);
        offGraphics.setColor(Color.green);
        offGraphics.drawRect(0,0,AppletWidth-1,AppletHeight-1);
        for(int i=0; i<NumberofStates; i++) {
          offGraphics.setColor(Color.blue);
          offGraphics.fillRect(xMargin+i*IncrementofX-5,AppletHeight-5-
                   HamilHeight[i],15,5);
          offGraphics.setColor(Color.yellow);
          offGraphics.fillRect(xMargin+i*IncrementofX,AppletHeight-
                     PulseHeight[i]-5-HamilHeight[i],5,PulseHeight[i]);
        }
        offGraphics.setColor(Color.orange);
        offGraphics.drawLine(xMargin+5,AppletHeight-5,xMargin+
                   IncrementofX,AppletHeight-HamilHeight[1]);
        offGraphics.drawLine(xMargin+(NumberofStates-1)*IncrementofX,
                   AppletHeight-5,
                   xMargin+(NumberofStates-2)*IncrementofX,
                   AppletHeight-HamilHeight[NumberofStates-2]);
        for(int i=1; i<NumberofStates-2; i++) {
          offGraphics.drawArc(xMargin+i*IncrementofX,AppletHeight-
                    HamilHeight[i]-20,IncrementofX,IncrementofX,10,160);
        }
        paint(onGraphics);
      }
      try { Thread.currentThread().sleep(delay); }
      catch (Exception e) {}
    }
    myThread=null;
  }

  public void update (Graphics g) { paint(g); }
  public void paint(Graphics g) {
    if(offImage!=null) {
      g.drawImage(offImage,0,0,AppletWidth,AppletHeight,null);
    }
  }
  

  private void getVibrations () {
    String parameter;
    double Amp = 2.0;
//    parameter = getParameter("period");
//    if (parameter!=null) { Period = Integer.parseInt(parameter); }
    parameter = getParameter("amplitude");
    if (parameter!=null) { Amp = Integer.parseInt(parameter); }
    for(int i=0; i<NumberofStates; i++) {
      if(i>0 && i<NumberofStates-1) {
        Amplitude[i] = Amp;
        Frequency[i] = BarrierHeight/Period;
      }
      else {
        Amplitude[i] = 0.0;
        Frequency[i] = 0.0;
      }
    }
  }

  private void getHamiltonian () {
    for(int i=0; i < NumberofStates; i++) {
      for(int j=0; j < NumberofStates; j++) {
        if(i==0 && j==0) {
          Hamiltonian[i][j] = 0.0;
        }
        if(i==NumberofStates-1 && j==NumberofStates-1) {
          Hamiltonian[i][j] = 0.0;
        }
        if(i==j && i!=0 && i!=NumberofStates-1) {
          Hamiltonian[i][j] = BarrierHeight;
        }
        if(Math.abs(i-j) > 0 && Math.abs(i-j) <= 1) {
          Hamiltonian[i][j] = 1.0;
        }
      }
    }
    if(printHamiltonian) {
      System.out.println("Hamiltonian Matrix");
      for (int i=0; i<NumberofStates; i++) {
        for (int j=0; j<NumberofStates; j++) {
          System.out.print(Hamiltonian[i][j]+"  ");
        }
        System.out.println();
      }
    }
  }

  private void resetPsi() {
    for(int i=0; i<NumberofStates; i++) {
      if(i==0) { 
        Psi[i] = new Complex(1,0); 
      }
      else { 
        Psi[i] = new Complex(0,0); 
      }
      if(printInitialPsi) {
        System.out.println(i + " " + Psi[i].Real() + " " + Psi[i].Imaginary());
      }
    }
  }

  private void RungeKutta() {

    Complex z,half,two,onesixth,z2,z3;
    double sum = 0.0;
    half = new Complex(0.5,0.0);
    two = new Complex(2.0,0.0);
    onesixth = new Complex(1.0/6.0,0.0);
    Complex[][] complexH = new Complex[NumberofStates][NumberofStates];
    Complex[] F1 = new Complex[NumberofStates];
    Complex[] F2 = new Complex[NumberofStates];
    Complex[] F3 = new Complex[NumberofStates];
    Complex[] F4 = new Complex[NumberofStates];
    Complex[] temp = new Complex[NumberofStates];
    
    for(int i=0; i<NumberofStates; i++) {
      for(int j=0; j<NumberofStates; j++) {
        if(i==j) { 
          complexH[i][j] = new Complex(0,-Steplength*
          (Hamiltonian[i][j]+Amplitude[i]*Math.sin(Frequency[i]*time)));
        }
        else {
          complexH[i][j] = new Complex(0,-Steplength*Hamiltonian[i][j]);
        }
      }
    }

// loop begins;

    IndexofStep ++;
    
    for(int i=0; i<NumberofStates; i++) {
      F1[i] = new Complex(0,0);
      for(int j=0; j<NumberofStates; j++) {
        z=complexH[i][j].Times(Psi[j]);
        F1[i]=F1[i].Plus(z);
      }
    }
    for(int i=0; i<NumberofStates; i++) {
      z=half.Times(F1[i]);
      temp[i]=Psi[i].Plus(z);
    }

    for(int i=0; i<NumberofStates; i++) {
      F2[i] = new Complex(0,0);
      for(int j=0; j<NumberofStates; j++) {
        z=complexH[i][j].Times(temp[j]);
        F2[i]=F2[i].Plus(z);
      }
    }
    for(int i=0; i<NumberofStates; i++) {
      z=half.Times(F2[i]);
      temp[i]=Psi[i].Plus(z);
    }

    for(int i=0; i<NumberofStates; i++) {
      F3[i] = new Complex(0,0);
      for(int j=0; j<NumberofStates; j++) {
        z=complexH[i][j].Times(temp[j]);
        F3[i]=F3[i].Plus(z);
      }
    }
    for(int i=0; i<NumberofStates; i++) {
      temp[i]=Psi[i].Plus(F3[i]);
    }

    for(int i=0; i<NumberofStates; i++) {
      F4[i] = new Complex(0,0);
      for(int j=0; j<NumberofStates; j++) {
        z=complexH[i][j].Times(temp[j]);
        F4[i]=F4[i].Plus(z);
      }
    }

    for(int i=0; i<NumberofStates; i++) {
      z2=two.Times(F2[i]);
      z3=two.Times(F3[i]);
      z=z2.Plus(z3);
      z2=z.Plus(F1[i]);
      z3=z2.Plus(F4[i]);
      z=onesixth.Times(z3);
      temp[i]=Psi[i].Plus(z);
      Psi[i]=temp[i];
    }

    time += Steplength;
    sum=0.0;
    for (int i=0; i<NumberofStates; i++) {
      Occupancy[i] = Psi[i].MagnitudeSquare();
      sum += Occupancy[i];
      PulseHeight[i] = (int) (Occupancy[i]*AppletWidth);
      HamilHeight[i] = (int) ((Hamiltonian[i][i]+
                       Amplitude[i]*Math.sin(Frequency[i]*time))*15);
    }
    if(IndexofStep%outputInterval==0) {    
      System.out.println(time + " " + Occupancy[0]
                              + " " + Occupancy[NumberofStates-1]
                              + " " + sum
                              + " " + Period);
    }

// loop ends;

  }

  private void Converter () {
    for (int i=0; i<NumberofStates; i++) {
      Occupancy[i] = Psi[i].MagnitudeSquare();
      PulseHeight[i] = (int) (Occupancy[i]*AppletWidth);
      HamilHeight[i] = (int) ((Hamiltonian[i][i]+
                       Amplitude[i]*Math.sin(Frequency[i]*time))*15);
    }
  }

  public boolean mouseDown(Event e, int x, int y) {
    if(threadSuspended) {
      myThread.resume();
    }
    else {
      myThread.suspend();
    }
    threadSuspended = !threadSuspended;
    return true;
  }

  public boolean action ( Event event, Object arg) {
    if(event.target instanceof Choice) {
      if(arg.equals("Highest")) {
        Period=0.5;
        getVibrations();
      }
      if(arg.equals("Higher")) {
        Period=1.0;
        getVibrations();
      }
      if(arg.equals("High")) {
        Period=10.0;
        getVibrations();
      }
      if(arg.equals("Middle")) {
        Period=50.0;
        getVibrations();
      }
      if(arg.equals("Low")) {
        Period=100.0;
        getVibrations();
      }
    }
    if(arg.equals("Reset")) { 
      time=0.0;
      IndexofStep=0;
      resetPsi();
      Converter();
      repaint();
      myThread.suspend(); 
    }
    if(arg.equals("Start")) {
      time=0.0;
      IndexofStep=0;
      delay=0;
      beginIt=true;
      myThread.resume();
    }
    return true;
  }

}

