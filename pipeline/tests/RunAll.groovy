import org.junit.internal.TextListener;
import org.junit.extensions.cpsuite.ClasspathSuite;
import org.junit.extensions.cpsuite.ClasspathSuite.ClassnameFilters;
import org.junit.runner.JUnitCore;
import org.junit.runner.RunWith;

@RunWith(ClasspathSuite.class)
@ClassnameFilters([".*Test"])
class RunAll {

    static void main(String[] args) {
        JUnitCore junit = new JUnitCore();
        junit.addListener(new TextListener(System.out));
        junit.run(RunAll.class);
    }

}
