//package specserve;

//init with $env:Path += ";C:\Program Files\Java\jdk1.8.0_171\bin;C:\Program Files\Ocean Optics\OmniDriver\OOI_HOME"
//	 javac -cp "C:\Users\optik\Documents\Code\Ocean\OmniDriver\OOI_HOME\OmniDriver.jar;
//    C:\Users\optik\Documents\Code\Ocean\OmniDriver\OOI_HOME\jSerialComm-1.3.11.jar" .\specserve2.java
//javac -cp "C:\Program Files\Ocean Optics\OmniDriver\OOI_HOME\OmniDriver.jar;C:
//\Program Files\Ocean Optics\OmniDriver\jSerialComm-2.0.2.jar" .\specserve2.java

import java.io.*;
import javax.xml.ws.*;
import javax.xml.ws.http.*;
import javax.xml.ws.handler.*;
import javax.xml.transform.*;
import javax.xml.transform.stream.*;
import javax.annotation.Resource;

import java.util.regex.*;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;
import java.util.Map;

import java.util.concurrent.TimeUnit;

import com.oceanoptics.omnidriver.api.wrapper.Wrapper;
import com.oceanoptics.omnidriver.features.gpio.GPIO;
import com.oceanoptics.omnidriver.features.detectortemperature.*;
import com.oceanoptics.omnidriver.features.thermoelectric.*;

import com.fazecast.jSerialComm.*;

@WebServiceProvider
@ServiceMode(value = Service.Mode.PAYLOAD)
public class specserve2 implements Provider<Source> {
    static int cnt;
    static Wrapper wrapper;
    static GPIO gpio;
	static int bin;
    Map<String,Integer> input;
    @Resource
    protected WebServiceContext wsctx;

    public static String data (int limit) {
        List<Integer> numbers = new ArrayList<>();
        for (int i=0;i<limit;i+=1) {
            numbers.add(5 + (int)(Math.random() * 11));
        }
        String out="";
        for (int j: numbers) {
            out=out.concat(String.format("<i>%d</i>",j));
        }
        return out;
    }


    public String expose (int intime,int chsel,int aver) {
      //int integrationTime = 1000;
        int nchan=wrapper.getWrapperExtensions().getNumberOfEnabledChannels(0);
        int maxPixels = 20000;
        int printSpec=20000;
        double[] wavelengths, spectralData;          // arrays of doubles to hold the spectral data
        int numberOfPixels;
        String rep="";
        double specSum=0;
        rep=rep+String.format("<exposure chan=\"%d\">%d</exposure> ",chsel,intime);
        for (int j = 0; j < nchan; j++){
            if ((chsel>=0) && (j!=chsel)) {continue;}
            // get a spectrum from the spectrometers:
            rep=rep+String.format("<data chan=\"%d\">",j);
            numberOfPixels = wrapper.getNumberOfPixels(0,j);          // gets the number of pixels in the first spectrometer.
            if ((maxPixels>0) && (numberOfPixels>maxPixels)) {numberOfPixels=maxPixels;}

            if (intime<=0) {
                wavelengths = wrapper.getWavelengths(0,j);                // gets the wavelengths of the first spectrometer and save them in a double array
                for (int i = 0; i < numberOfPixels; i++){
                    rep=rep+String.format("%5.3f ", wavelengths[i]);
                }

            }  else {
                wrapper.setIntegrationTime(0,j, intime*1000);         // sets the integration time of the first spectrometer
                if (aver>0) {wrapper.setScansToAverage(0,j, aver);}
                spectralData = wrapper.getSpectrum(0,j);                  // gets the spectrum from the first spectrometer and saves it to a double array
                System.out.printf("# channel "+j+"  : \n");
                // loop for printing the spectral data to the screen:
                specSum=0;
				double binSum=0;
                for (int i = 0; i < numberOfPixels; i++){
                    if (i<printSpec)  {
						if (this.bin>0) {
							binSum+=spectralData[i];
							if ((i+1)%this.bin==0) {
								rep=rep+String.format("%5.3f ", binSum/this.bin);
								binSum=0;
							}
						} else {
							rep=rep+String.format("%5.3f ", spectralData[i]);
						}
					}
                    specSum+=spectralData[i];
                }
            }
            rep=rep+String.format("</data><sum chan=\"%d\">%.0f</sum>",j,specSum);
        }
        return rep;
    }

	public String exposemul (int []intime, int aver) {
    String rep="";
    int j=0;
		for (int i=0; i<intime.length; i++) {
			if (intime[i]>0) {
        rep=rep+expose(intime[i],i,aver);
        j+=1;
			}
    }
    rep=rep+String.format("<p>multichanel %d</p>",j);
		return rep;
  }

  public String expoline(int nstep,int dstep,int intime,int chsel,int aver,char dir) //throws InterruptedException
  {
    String rep="";
    SerialPort comPort = SerialPort.getCommPorts()[0];
    //SerialPort comPort = SerialPort.getCommPort("COM4");
    comPort.setBaudRate(250000);
    comPort.setComPortTimeouts(1,1,1);
    comPort.openPort();
    int pos=20;
    //byte[] wBuffer;
    String ostr= new String("M114\n"); //get position
    try {
      comPort.writeBytes(ostr.getBytes(),ostr.length());
      TimeUnit.SECONDS.sleep(1);
      // ma vratit ok C: X:0.00 Y:0.00 Z:0.00 E:0.00
    } catch (InterruptedException e) {};
	int nbytes=100;//comPort.bytesAvailable();
	byte[] readBuffer;
	int numRead;
	if (nbytes>0) {
        readBuffer = new byte[nbytes];
		numRead = comPort.readBytes(readBuffer, readBuffer.length);
		if (numRead>1) {
		  String inp1= new String(readBuffer);
		  String patt=dir+":\\d+";
		  Pattern pattern = Pattern.compile(patt,Pattern.CASE_INSENSITIVE);
		  Matcher matcher = pattern.matcher(inp1);
		  if (matcher.find()) {
			try {
			  pos=Integer.parseInt(inp1.substring(matcher.start()+2,matcher.end()));
			  System.out.println(inp1+": position "+pos);
			} catch (NumberFormatException e) {
			  System.out.println(": problem parsing input "+inp1);
			}
		  } else {
			System.out.println(inp1+":parse failed");
		  }
		}
	}
    int i;
    for (i=0;i<nstep;i+=1) {
      pos+=dstep;
      //ostr= String.format("G1 Y%i\n",pos);
      ostr= new String("G1 "+dir);
      ostr=ostr+pos+"\n";
      comPort.writeBytes(ostr.getBytes(),ostr.length());
      rep=rep+expose(intime,chsel,aver);
      rep+="<pos>"+Integer.toString(pos)+"</pos>";
    }

    try {
    if (1==0) {
       while (true)
       {
          while (comPort.bytesAvailable() == 0)
            Thread.sleep(20);

          readBuffer = new byte[comPort.bytesAvailable()];
          numRead = comPort.readBytes(readBuffer, readBuffer.length);
          System.out.println("Read " + numRead + " bytes.");
       }
    }
    } catch (Exception e) { e.printStackTrace(); }
    comPort.closePort();
    return rep;
  }

	public Source invoke(Source request) {
        MessageContext ctx=wsctx.getMessageContext();
        String query=(String) ctx.get(MessageContext.QUERY_STRING);
        int eqpos;
        String rep=String.format("<p>Hello there [%d]!",cnt);
        Source reply;
        input = new HashMap<String,Integer>();
        char dir='Y';
        if (query!=null) {
          for (String asgval: query.split("&")) {
              eqpos=asgval.indexOf("=");
              if (eqpos>=0) {
                  System.out.println("found "+asgval.substring(0,eqpos));
                  if (asgval.substring(0,eqpos).equals("dstep")) {
                    dir=asgval.charAt(eqpos+1);
                    input.put(asgval.substring(0,eqpos),Integer.parseInt(asgval.substring(eqpos+2)));
                  } else {
                    input.put(asgval.substring(0,eqpos),Integer.parseInt(asgval.substring(eqpos+1)));
                  }
              }
          }
        } else {
          rep=rep+" check the usage ";
        }
        cnt=cnt+1;
        int intime=10;
        int aver=10;
        int specsel=0;
        int chsel=-1;

        int nstep=1;
        int dstep=10;
        //int ipos=0;
        int shutterbit=3;
        //assert  shutterbit<gpio.getTotalBits()
        int shutterval=0;
        String key5="shut";
        if (input.containsKey(key5)) {
            try {
                shutterval=input.get(key5);
                if (shutterval>=1) {
                  //System.out.println("switching on");
                  wrapper.setStrobeEnable(specsel,1);}
                else {
                  //System.out.println("switching off");
                  wrapper.setStrobeEnable(specsel,0);
                }
               //gpio.setValueBit(shutterbit,shutterval==1); //nebo setMuxBit
            } catch (Exception ee) {}
        }
        key5="setbit";
        if (input.containsKey(key5)) {
            try {
               shutterbit=input.get(key5);
               System.out.println("setting bit "+shutterbit);
               gpio.setMuxBit(shutterbit,true); //nebo setMuxBit
            } catch (Exception ee) {}
        }
        key5="delbit";
        if (input.containsKey(key5)) {
            shutterbit=input.get(key5);
            //gpioBitSet = gpio.getValueBits();
            //shutterval= gpioBitSet.get(shutterbit);
            System.out.println("bit "+shutterbit+" is "+shutterval);
            try {
               gpio.setMuxBit(shutterbit,false); //nebo setMuxBit
            } catch (Exception ee) {}
        }
        String key6="temp";
        double btemp=0;
        try {
            if (input.containsKey(key6)) {
               if (wrapper.isFeatureSupportedThermoElectric(specsel)) {
                  btemp=input.get(key6);
                  ThermoElectricWrapper tdet=wrapper.getFeatureControllerThermoElectric(specsel);
        		      //DetectorTemperature det=wrapper.getFeatureControllerDetectorTemperature(specsel);
                  if (tdet!=null) {
                    tdet.setDetectorSetPointCelsius(btemp);
                    btemp=tdet.getDetectorTemperatureCelsius();
                    System.out.println("# temperature "+btemp);
                  } else {
                    System.out.println("# temperature unknown");
                  }
      		       }
               }
             } catch (Exception e) { e.printStackTrace(); }
         String key7="corr";
         if (input.containsKey(key7)) {
              int nlin=1;
              int icors;
              icors=input.get(key7);
              System.out.println("# apply corrections");
              wrapper.setCorrectForDetectorNonlinearity(specsel, icors%2);
              wrapper.setCorrectForElectricalDark(specsel, icors/2);
              rep=rep+" correcting...";
         }
         String key8="step";
         if (input.containsKey(key8)) {
           nstep=input.get(key8);
           key8="dstep";
           dstep=input.get(key8);
         }
         String key9="bin";
         if (input.containsKey(key9)) {
           this.bin=input.get(key9);
         }
          String key0="avg";
          if (input.containsKey(key0)) {
               aver=input.get(key0);
          }
          String key1="exp";
          if (input.containsKey(key1)) {
               intime=input.get(key1);
               rep=rep+" measuring...";
          }
          String key2="chan";
          if (input.containsKey(key2)) {
               chsel=input.get(key2);
          }
		key1="exp0";
    if (chsel<0 && input.containsKey(key1)) {//multichannel measurement
      int[] inmult;
			inmult=new int[10];
      int i;
      int j=0;
			for (i=0;i<10;i+=1) {
        inmult[i]=0;
				//key1=String.format("exp%i",i);
				key1="exp"+Integer.toString(i);
				if (input.containsKey(key1)) {
            inmult[i]=input.get(key1);
            j+=1;
				} else {
					//inmult[i]=0;
					break;
				}
			}
			//if (j>0) {
        rep="<p>Hello multichannels!";
        reply = new StreamSource(new StringReader(rep + exposemul(Arrays.copyOfRange(inmult,0,i),aver) + "</p>"));
        return reply;
      //}
		}
    if (nstep>1) {
        System.out.printf("# scanning "+nstep+" pts \n");
        reply = new StreamSource(new StringReader(rep + expoline(nstep,dstep,intime,chsel,aver,dir) + "</p>"));
      }
    else
        reply = new StreamSource(new StringReader(rep + expose(intime,chsel,aver) + "</p>"));
    if (shutterval==2) wrapper.setStrobeEnable(specsel,0);
    return reply;
  }


  public static void main(String[] args) throws InterruptedException {

        String address = "http://127.0.0.1:";
        if (args.length >0) {address+=args[0]+"/";}
        else {address+="8080/";}
        Endpoint.create(HTTPBinding.HTTP_BINDING, new specserve2()).publish(address);

        wrapper = new Wrapper();             // Create a wrapper object

        int numberOfSpectrometers;   // variables for the program
        int integrationTime = 1000;                  // set the integration time
        String serialNumber;                         // variable to get the serial number
                                                     // and wavelengths.
        System.out.println("Starting");
        numberOfSpectrometers = wrapper.openAllSpectrometers(); // opens all of the spectrometers and returns the number of spectrometers found
        if (numberOfSpectrometers < 1){              // Check for any spectrometers
            System.out.println("There are no spectrometers attached to the computer");
            return;
        } else {
            System.out.println("Something found");
        }
        serialNumber = wrapper.getSerialNumber(0);              // gets the serial number from the first spectrometer
        System.out.println("Serial Number: " + serialNumber);   // prints the serial number to the screen
        gpio = wrapper.getFeatureControllerGPIO(0);
        int nchan=wrapper.getWrapperExtensions().getNumberOfEnabledChannels(0);
        System.out.println("found "+nchan+" channels");

        System.out.println("Service running at " + address);
        System.out.println("Type [CTRL]+[C] to quit!");
        //System.out.println(data(10));
        Thread.sleep(Long.MAX_VALUE);

		//gpio.setMuxBit(7,True)
		//gpio.setValueBit(7,True)

    }
}
