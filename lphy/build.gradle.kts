import java.text.SimpleDateFormat
import java.util.Calendar

plugins {
    `java-library`
}

group = "lphy"
version = "1.1-SNAPSHOT"


dependencies {
    // required in test
    api("org.antlr:antlr4-runtime:4.8")
    api("org.apache.commons:commons-math3:3.6.1")
    api("org.apache.commons:commons-lang3:3.10")
    // in maven
    api("org.scilab.forge:jlatexmath:1.0.7")
    api("org.scilab.forge:jlatexmath-font-greek:1.0.7")
    api("org.scilab.forge:jlatexmath-font-cyrillic:1.0.7")
    api("net.steppschuh.markdowngenerator:markdowngenerator:1.3.1.1")

    testImplementation("junit:junit:4.13.2")
//    testRuntimeOnly("org.junit.jupiter:junit-jupiter-engine:4.13")

    // not in maven
    api(files("libs/jebl-3.0.1.jar"))
    //implementation(fileTree("lib") { exclude("junit-*.jar") })

}

java {
    sourceCompatibility = JavaVersion.VERSION_16
    targetCompatibility = JavaVersion.VERSION_16
}

// overwrite compileJava to use module-path
// overwrite compileJava to pass dependencies to tests
tasks.compileJava {
    // use the project's version or define one directly
    options.javaModuleVersion.set(provider { project.version as String })

    println("Java version used is ${JavaVersion.current()}.")

    doFirst {
        println("CLASSPATH IS ${classpath.asPath}")
        options.compilerArgs = listOf("--module-path", classpath.asPath)
        classpath = files()
    }
}

var calendar:Calendar? = Calendar.getInstance()
var formatter = SimpleDateFormat("dd-MMM-yyyy HH:mm:ss")

tasks.jar {
    manifest {
        // shared attr in the root build
        attributes(
            "Implementation-Title" to "LPhy",
            "Implementation-Version" to archiveVersion,
            "Built-Date" to formatter.format(calendar?.time)
        )
    }
}

tasks.test {
    useJUnit()
    // useJUnitPlatform()
    maxHeapSize = "1G"
}
