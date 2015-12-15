import barrypitman.junitXmlFormatter.*
import org.junit.extensions.cpsuite.ClasspathSuite
import org.junit.extensions.cpsuite.ClasspathSuite.ClassnameFilters;;
import org.junit.internal.TextListener;
import org.junit.runner.JUnitCore;
import org.junit.runner.RunWith;

@RunWith(ClasspathSuite.class)
@ClassnameFilters([".*Test"])
class RunAll {

    static void main(String... args) {
        AntXmlRunListener runListener = new AntXmlRunListener();
        try {
            runListener.setOutputStream(new FileOutputStream(new File("TEST-Result.xml")));
        }catch (FileNotFoundException msg){
            System.err.println("Test result report cannot be generated.");
        }
        JUnitCore junit = new JUnitCore();
        junit.addListener(runListener);
        junit.addListener(new TextListener(System.out));
        junit.run(RunAll.class);
    }

}
