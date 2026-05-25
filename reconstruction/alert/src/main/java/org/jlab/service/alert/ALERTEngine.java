package org.jlab.service.alert;

import ai.djl.translate.TranslateException;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.HashMap;

import org.jlab.clas.reco.ReconstructionEngine;
import org.jlab.clas.swimtools.Swim;
import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.geom.detector.alert.ATOF.AlertTOFFactory;
import org.jlab.geom.detector.alert.ATOF.AlertTOFDetector;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.rec.alert.TrackMatchingAI.ModelTrackMatching;
import org.jlab.rec.alert.AIPID.ModelPrePID;
import org.jlab.rec.alert.banks.RecoBankWriter;
import org.jlab.rec.alert.projections.TrackProjector;
import org.jlab.rec.atof.hit.ATOFHit;
import org.jlab.rec.ahdc.KalmanFilter.RadialKFHit;
import org.jlab.rec.ahdc.KalmanFilter.KalmanFilter;
import org.jlab.rec.ahdc.Hit.Hit;
import org.jlab.geom.detector.alert.AHDC.AlertDCDetector;
import org.jlab.geom.detector.alert.AHDC.AlertDCFactory;
import org.jlab.geom.detector.alert.AHDC.AlertDCWireIdentifier;
import org.jlab.rec.ahdc.Track.Track;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.jlab.clas.pdg.PDGDatabase;
import org.jlab.clas.pdg.PDGParticle;
import java.util.logging.Logger;
import org.jlab.rec.alert.constants.CalibrationConstantsLoader;



import ai.djl.util.Pair;
import org.jlab.rec.alert.AIPID.PrePIDResult;


/** 
 * <h1>ALERTEngine reconstruction service.</h1>
 *
 * @author  Whit Armstrong
 * @author  Noemie Pilleux
 * @since   2025-04-03
 */
public class ALERTEngine extends ReconstructionEngine {

    /**
     * ALERT Engine output bank writer.
     * 
     * @see RecoBankWriter
     *
     * <h3>Output banks</h3>
     * <ul>
     * <li> Track Projection @see TrackProjector</li>
     * </ul>
     *
     */
    private RecoBankWriter rbc;
    static final Logger LOGGER = Logger.getLogger(ModelPrePID.class.getName());
    AlertTOFDetector ATOF; // ALERT ATOF detector
    private AlertDCDetector AHDC; // ALERT AHDC detector

    /**
     *  Current run number being processed.
     *  TODO: why atomic here and nowhere else? 
     */
    private final AtomicInteger run = new AtomicInteger(0);

    private double b; //Magnetic field

    private ModelTrackMatching modelTrackMatching;
    private ModelPrePID modelPrePID;

    public void setB(double B) {
        this.b = B;
    }
    public double getB() {
        return b;
    }

    /**
     * ALERTEngine service c'tor. 
     */
    public ALERTEngine() {
        super("ALERT", "whit,ouillon,pilleux", "0.1");
    }

    /** 
     * ALERTEngine initialization.
     * Creates the RecoBankWriter and checks for various yaml flags.
     * TODO: document flags
     */
    @Override
    public boolean init() {

        rbc = new RecoBankWriter();

        modelTrackMatching = new ModelTrackMatching();
        modelPrePID = new ModelPrePID();

        AlertTOFFactory factory = new AlertTOFFactory();
        DatabaseConstantProvider cp = new DatabaseConstantProvider(11, "default");
        ATOF = factory.createDetectorCLAS(cp);

        /// --- apply current layer alignment
        double[] layer_angles_start = {0.8609, 1.0181, 0.5654, 0.7998, 0.3913, 0.5151, 0.2749, 0.5057};
        double[] layer_angles_end   = {0.8412, 0.8157, 0.4084, 0.7939, 0.4747, 0.5086, 0.3518, 0.2534};
        double[] wire_angles_start = layerAngles2WireAngles(layer_angles_start);
        double[] wire_angles_end   = layerAngles2WireAngles(layer_angles_end);
        AHDC = generateAhdcGeometry(wire_angles_start, wire_angles_end);
        /// --- end alignment

        //AHDC = (new AlertDCFactory()).createDetectorCLAS(new DatabaseConstantProvider());

        if(this.getEngineConfigString("Mode")!=null) {
            //if (Objects.equals(this.getEngineConfigString("Mode"), Mode.AI_Track_Finding.name()))
            //    mode = Mode.AI_Track_Finding;
        }

        // Requires calibration constants
        String[] alertTables = new String[] {
            "/calibration/alert/ahdc/time_offsets",
            "/calibration/alert/ahdc/time_to_distance",
            "/calibration/alert/ahdc/raw_hit_cuts",
            "/calibration/alert/atof/effective_velocity",
            "/calibration/alert/atof/time_walk",
            "/calibration/alert/atof/attenuation",
            "/calibration/alert/atof/time_offsets",
            "/calibration/alert/ahdc/gains",
            "/calibration/alert/ahdc/time_over_threshold"
        };
        requireConstants(Arrays.asList(alertTables));
        return true;
    }

    /**
     * Process Event.
     * Main method called to process event data.
     *
     * <ul>
     * <li> Check for AHDC and ATOF banks </li>
     * <li> Project track to ATOF</li>
     * </ul>
     */
    @Override
    public boolean processDataEvent(DataEvent event) {
        return processDataEvent(event, AHDC);
    }

    long computing_time = 0;

    /** Retrun computing time in nanoseconds */
    public long getComputingTime() { return computing_time;}

    public boolean processDataEvent(DataEvent event, AlertDCDetector AHDCdet) {

        long start_time = System.nanoTime();

        if (!event.hasBank("AHDC::adc")) 
            return false;
        // if (!event.hasBank("ATOF::tdc")) 
        //     return false;

        if (!event.hasBank("RUN::config")) {
            return true;
        }

        DataBank runBank = event.getBank("RUN::config");

        int newRun = runBank.getInt("run", 0);
        if (newRun == 0) {
            return true;
        }

        if (run.get() == 0 || (run.get() != 0 && run.get() != newRun)) {
            run.set(newRun);
        }
      
        ///////////////////////////////////////////
        /// Kalmam Filter
        /// ///////////////////////////////////////
        
        /// Pre conditions
        if (!event.hasBank("AHDC::track")) {return false;}
        if (!event.hasBank("AHDC::hits")) {return false;}

        /// tmp: misalignement with respect to the center of the AHDC (mm)
        /// the AHDC is at the center of the solenoid
        // double clas_alignement = +30;
        // double atof_alignement = -32.7;
        

        /// Read the electron vertex
        double vz_electron = 0;
        double[] vz_error2 = {0.09, 1e10, 1e10}; // mm^2, radians^2, mm^2, error on r, phi, z
        boolean IsVertexDefined = false;
        if (event.hasBank("REC::Particle")) {
            DataBank recBank = event.getBank("REC::Particle");
            for (int row = 0; row < recBank.rows(); row++) {
                if (recBank.getInt("pid", row) == 11) {
                    vz_electron = 10*recBank.getFloat("vz",row); // conversion in mm
                    IsVertexDefined = true;
                    

                    //double px = recBank.getFloat("px",row);
                    //double py = recBank.getFloat("py",row);
                    //double pz = recBank.getFloat("pz",row);
                    //double p = Math.sqrt(px*px+py*py+pz*pz);
                    //double theta = Math.acos(pz/p);
                    
                    // set the resolutions on r and z! to be done
                    vz_error2[0] = 0.09; // should depend on p and theta
                    vz_error2[2] = 64; // should depend on p and theta

                    break; // only look at the first electron
                }
            }
        }

        /// Read the list of tracks/hits from the banks AHDC::track and AHDC::hits
        DataBank trackBank = event.getBank("AHDC::track");
        DataBank hitBank = event.getBank("AHDC::hits");
        ArrayList<Track> AHDC_tracks = new ArrayList<>();
        for (int row = 0; row < trackBank.rows(); row++) {
            int trackid = trackBank.getInt("trackid", row);
            ArrayList<Hit> AHDC_hits = new ArrayList<>();
            for (int hit_row = 0; hit_row < hitBank.rows(); hit_row++) {
                if(trackid == hitBank.getInt("trackid", hit_row)) {
                    int id = hitBank.getShort("id", hit_row);
                    int superlayer = hitBank.getByte("superlayer", hit_row);
                    int layer = hitBank.getByte("layer", hit_row);
                    int wire = hitBank.getInt("wire", hit_row);
                    int adc = hitBank.getInt("adc", hit_row);
                    double doca = hitBank.getDouble("doca", hit_row);
                    double time = hitBank.getDouble("time", hit_row);
                    double tot = hitBank.getDouble("timeOverThreshold", hit_row);
                    // warning : adc is the calibrated one, we need the adc for the Kalman filter !
                    int sector = 1; // constant value
                    int key_value = sector*10000 + (superlayer*10 + layer)*100 + wire;
                    double raw_adc = adc;
                    double[] adc_gain = CalibrationConstantsLoader.AHDC_ADC_GAINS.get(key_value);
                    if (adc_gain != null) raw_adc = adc/adc_gain[0];
                    //System.out.println("adc : " + adc + " raw_adc : " + raw_adc);
                    Hit hit = new Hit(id, superlayer, layer, wire, doca, raw_adc, time);
                    hit.setWirePosition(AHDCdet);
                    hit.setTrackId(trackid);
                    hit.setADC(adc);
                    hit.setToT(tot);
                    AHDC_hits.add(hit);
                }
            }
            if (AHDC_hits.size() > 0) {
                Collections.sort(AHDC_hits); // sorted following the compareTo() method in Hit.java
                Track track = new Track(AHDC_hits);
                // Initialise the position and the momentum using the information of the AHDC::track
                // position : mm
                // momentum : MeV
                double x = trackBank.getFloat("x", row);
                double y = trackBank.getFloat("y", row);
                double z = trackBank.getFloat("z", row);
                double px = trackBank.getFloat("px", row);
                double py = trackBank.getFloat("py", row);
                double pz = trackBank.getFloat("pz", row);
                double[] vec = {x, y, z, px, py, pz};
                track.setPositionAndMomentumVec(vec);
                track.set_trackId(trackid);
                AHDC_tracks.add(track);
            }
        }

        // /// Associate the electron vertex (the beamline hit) to each track
        boolean IsMC = event.hasBank("MC::Particle");
        double vz_constraint = vz_electron + (IsMC ? 0 : clas_alignement); // we don't have the misalignment in simulation
        for (Track track : AHDC_tracks) {
            RadialKFHit hit_beam = new RadialKFHit(0, 0, vz_constraint);
            RealMatrix measurementNoise = new Array2DRowRealMatrix(
											new double[][]{
												{vz_error2[0], 0.0000      , 0.0000},
												{0.0000      , vz_error2[1], 0.0000},
												{0.0000      , 0.0000      , vz_error2[2]}
											});//3x3;
			hit_beam.setMeasurementNoise(measurementNoise);
            track.setBeamlineHit(hit_beam);
        }

        /// Intialise the Kalman Filter
        double magfieldfactor = runBank.getFloat("solenoid", 0);
        double magfield = 50*magfieldfactor;
        PDGParticle proton = PDGDatabase.getParticleById(2212);
        KalmanFilter KF = new KalmanFilter(proton, Niter);
        
        // ------tmp propagation without a fit
        //KF.propagationWithoutCorrection(AHDC_tracks, magfield, IsMC, event);
        // ----------
        
        KF.set_ATOF_detector(null); // Reference the ATOF geometry in the Kalman Filter
        KF.set_atof_alignement(atof_alignement);
        KF.set_vz_constraint(vz_constraint);
        KF.set_vertex_flag(IsVertexDefined);
        KF.setStepSize(stepper_size);

        /// Do a first propagation
        KF.propagation(AHDC_tracks, magfield, IsMC);

        /// Clean AHDC bad hits
        double sigma = 0.5; // mm
        for (Track track : AHDC_tracks) {
            ArrayList<Hit> AHDC_hits = track.getHits();
            Iterator<Hit> it = AHDC_hits.iterator();
            while (it.hasNext()) {
                Hit hit = it.next();
                if (Math.abs(hit.getResidual()) > 3*sigma) {
                    it.remove();
                }
            }
        }

        // /// Second propagation : each AHDC_tracks will be fitted
        KF.set_Niter(10);
        KF.propagation(AHDC_tracks, magfield, IsMC);

        /// write the AHDC::kftrack bank in the event
        event.removeBank("AHDC::kftrack");
        org.jlab.rec.ahdc.Banks.RecoBankWriter ahdc_writer = new org.jlab.rec.ahdc.Banks.RecoBankWriter();
        DataBank recoKFTracksBank   = ahdc_writer.fillAHDCKFTrackBank(event, AHDC_tracks);
        event.appendBank(recoKFTracksBank);
        // update the AHDC::hits bank : fill the residuals
        event.removeBank("AHDC::hits");
        ArrayList<Hit> AHDC_hits = new ArrayList<>();
        for (Track track : AHDC_tracks) {
            AHDC_hits.addAll(track.getHits());
        }     
        DataBank recoKFHitsBank = ahdc_writer.fillAHDCHitsBank(event, AHDC_hits);
        event.appendBank(recoKFHitsBank); // remark: only  hits associated to a track are saved

        long end_time = System.nanoTime();
        computing_time = end_time - start_time;
 

        return true;
    }

    double clas_alignement = +75;
    double atof_alignement = 0;

    public void set_atof_alignement(double _shift) {this.atof_alignement = _shift;}
	public void set_clas_alignement(double _shift) {this.clas_alignement = _shift;}

    private double stepper_size = 0.5;
	public void setStepSize(double _size) { stepper_size = _size;}
	public double getStepSize() { return stepper_size;}

    int Niter = 30;
    public void set_KF_Niter(int Niter) {this.Niter = Niter;}
	public int  get_KF_Niter() {return this.Niter;}

    public boolean processDataEventProjOnly(DataEvent event, AlertDCDetector AHDCdet) {

        if (!event.hasBank("AHDC::adc")) 
            return false;
        if (!event.hasBank("ATOF::tdc")) 
            return false;

        if (!event.hasBank("RUN::config")) {
            return true;
        }

        DataBank runBank = event.getBank("RUN::config");

        int newRun = runBank.getInt("run", 0);
        if (newRun == 0) {
            return true;
        }

        if (run.get() == 0 || (run.get() != 0 && run.get() != newRun)) {
            run.set(newRun);
        }

        
        // //Do we need to read the event vx,vy,vz?
        // //If not, this part can be moved in the initialization of the engine.
        // double eventVx=0,eventVy=0,eventVz=0; //They should be in CM
        // //Track Projector Initialisation with b field
        // Swim swim = new Swim();
        // float magField[] = new float[3];
        // swim.BfieldLab(eventVx, eventVy, eventVz, magField); 
        // this.b = Math.sqrt(Math.pow(magField[0],2) + Math.pow(magField[1],2) + Math.pow(magField[2],2));
        

        // TrackProjector projector = new TrackProjector();
        // projector.setB(this.b);
        // projector.projectTracks(event);
        // rbc.appendMatchBanks(event, projector.getProjections());

        ///////////////////////////////////////////
        /// Kalmam Filter
        /// ///////////////////////////////////////
        
        /// Read the list of tracks/hits from the banks AHDC::track and AHDC::hits
        DataBank trackBank = event.getBank("AHDC::track");
        DataBank hitBank = event.getBank("AHDC::hits");
        ArrayList<Track> AHDC_tracks = new ArrayList<>();
        for (int row = 0; row < trackBank.rows(); row++) {
            int trackid = trackBank.getInt("trackid", row);
            ArrayList<Hit> AHDC_hits = new ArrayList<>();
            for (int hit_row = 0; hit_row < hitBank.rows(); hit_row++) {
                if(trackid == hitBank.getInt("trackid", hit_row)) {
                    int id = hitBank.getShort("id", hit_row);
                    int superlayer = hitBank.getByte("superlayer", hit_row);
                    int layer = hitBank.getByte("layer", hit_row);
                    int wire = hitBank.getInt("wire", hit_row);
                    int adc = hitBank.getInt("adc", hit_row);
                    double doca = hitBank.getDouble("doca", hit_row);
                    double time = hitBank.getDouble("time", hit_row);
                    double tot = hitBank.getDouble("timeOverThreshold", hit_row);
                    // warning : adc is the calibrated one, we need the adc for the Kalman filter !
                    int sector = 1; // constant value
                    int key_value = sector*10000 + (superlayer*10 + layer)*100 + wire;
                    double raw_adc = adc;
                    double[] adc_gain = CalibrationConstantsLoader.AHDC_ADC_GAINS.get(key_value);
                    if (adc_gain != null) raw_adc = adc/adc_gain[0];
                    //System.out.println("adc : " + adc + " raw_adc : " + raw_adc);
                    Hit hit = new Hit(id, superlayer, layer, wire, doca, raw_adc, time);
                    hit.setWirePosition(AHDCdet); // key point
                    hit.setTrackId(trackid);
                    hit.setADC(adc);
                    hit.setToT(tot);
                    AHDC_hits.add(hit);
                }
            }
            if (AHDC_hits.size() > 0) {
                Collections.sort(AHDC_hits); // sorted following the compareTo() method in Hit.java
                Track track = new Track(AHDC_hits);
                // Initialise the position and the momentum using the information of the AHDC::track
                // position : mm
                // momentum : MeV
                double x = trackBank.getFloat("x", row);
                double y = trackBank.getFloat("y", row);
                double z = trackBank.getFloat("z", row);
                double px = trackBank.getFloat("px", row);
                double py = trackBank.getFloat("py", row);
                double pz = trackBank.getFloat("pz", row);
                double[] vec = {x, y, z, px, py, pz};
                track.setPositionAndMomentumVec(vec);
                track.set_trackId(trackid);
                AHDC_tracks.add(track);
            }
        }

        // /// Associate the electron vertex (the beamline hit) to each track
        boolean IsMC = event.hasBank("MC::Particle");
        

        /// Intialise the Kalman Filter
        double magfieldfactor = runBank.getFloat("solenoid", 0);
        double magfield = 50*magfieldfactor;
        PDGParticle proton = PDGDatabase.getParticleById(2212);
        
        KalmanFilter KF = new KalmanFilter(proton, Niter);
        
        // ------tmp propagation without a fit
        KF.set_clas_alignement(this.clas_alignement); //mm
        KF.set_atof_alignement(this.atof_alignement); // mm
        KF.propagationWithoutCorrection(AHDC_tracks, magfield, IsMC, event);
        KF.setStepSize(stepper_size);
        // ----------
        

        /// write the AHDC::kftrack bank in the event
        event.removeBank("AHDC::kftrack");
        org.jlab.rec.ahdc.Banks.RecoBankWriter ahdc_writer = new org.jlab.rec.ahdc.Banks.RecoBankWriter();
        DataBank recoKFTracksBank   = ahdc_writer.fillAHDCKFTrackBank(event, AHDC_tracks);
        event.appendBank(recoKFTracksBank);
        // update the AHDC::hits bank : fill the residuals
        event.removeBank("AHDC::hits");
        ArrayList<Hit> AHDC_hits = new ArrayList<>();
        for (Track track : AHDC_tracks) {
            AHDC_hits.addAll(track.getHits());
        }     
        DataBank recoKFHitsBank = ahdc_writer.fillAHDCHitsBank(event, AHDC_hits);
        event.appendBank(recoKFHitsBank); // remark: only  hits associated to a track are saved
 

        return true;
    }

    /**
     * ALERTEngine main.
     * TODO: needs good test.
     */
    public static void main(String[] args) {

        double starttime = System.nanoTime();

        int    nEvent     = 0;
        int    maxEvent   = 1000;
        int    myEvent    = 3;
        String inputFile  = "alert_out_update.hipo";
        String outputFile = "output.hipo";

        if (new File(outputFile).delete()) System.out.println("output.hipo is delete.");

        System.err.println(" \n[PROCESSING FILE] : " + inputFile);

        ALERTEngine en = new ALERTEngine();

        HipoDataSource reader = new HipoDataSource();
        HipoDataSync   writer = new HipoDataSync();

        en.init();

        reader.open(inputFile);
        writer.open(outputFile);

        while (reader.hasEvent() && nEvent < maxEvent) {
            nEvent++;
            // if (nEvent % 100 == 0) System.out.println("nEvent = " + nEvent);
            DataEvent event = reader.getNextEvent();

            // if (nEvent != myEvent) continue;
            // System.out.println("***********  NEXT EVENT ************");
            // event.show();

            en.processDataEvent(event);
            writer.writeEvent(event);

        }
        writer.close();

        System.out.println("finished " + (System.nanoTime() - starttime) * Math.pow(10, -9));
    }

    // to fix dependency issue
    void doNothing() {
        //Do we need to read the event vx,vy,vz?
        //If not, this part can be moved in the initialization of the engine.
        double eventVx=0,eventVy=0,eventVz=0; //They should be in CM
        //Track Projector Initialisation with b field
        Swim swim = new Swim();
        float magField[] = new float[3];
        swim.BfieldLab(eventVx, eventVy, eventVz, magField); 
        this.b = Math.sqrt(Math.pow(magField[0],2) + Math.pow(magField[1],2) + Math.pow(magField[2],2));
    }


    /**
     * Code copied from amon/java-utils/.../AhdcAlignmentAnalyser
     * 
     * Convert a layer by layer results to a wire by wire results. Idea: all the wires belonging to the same layer have the same modification.
     * @param layer_angles
     * @return a vector of 576 double containing the wire values
     */
    static public double[] layerAngles2WireAngles(double[] layer_angles) {
        double[] wire_angles = new double[576];
        for (int i = 0; i < 576; i++) {
            AlertDCWireIdentifier identifier = new AlertDCWireIdentifier(i);
            int num = AlertDCWireIdentifier.layer2number(identifier.getLayerId())-1;
            wire_angles[i] = layer_angles[num];
        }
        return wire_angles;
    }

    /**
     * Code copied from amon/java-utils/.../AhdcAlignmentAnalyser
     * Generate the AHDC geometry with specific correction angles
     * @param _wire_angles_start rotation to be applied to the start of the AHDC wires
     * @param _wire_angles_end rotation to be applied to the end of the AHDC wires
     * @return AHDC geometry
     */
    static AlertDCDetector generateAhdcGeometry(double[] _wire_angles_start, double[] _wire_angles_end) {
        AlertDCFactory factory = new AlertDCFactory();
        factory.setWireCorrectionAngles(_wire_angles_start, _wire_angles_end);
        return factory.createDetectorCLAS(new DatabaseConstantProvider());
    }
}
