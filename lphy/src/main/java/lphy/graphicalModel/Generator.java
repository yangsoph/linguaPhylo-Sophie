package lphy.graphicalModel;

import lphy.core.narrative.Narrative;
import lphy.parser.functions.ExpressionNode;
import net.steppschuh.markdowngenerator.link.Link;
import net.steppschuh.markdowngenerator.list.UnorderedList;
import net.steppschuh.markdowngenerator.text.Text;
import net.steppschuh.markdowngenerator.text.emphasis.BoldText;
import net.steppschuh.markdowngenerator.text.heading.Heading;

import java.lang.annotation.Annotation;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Type;
import java.util.*;

/**
 * A generator generates values, either deterministically (DeterministicFunction) or stochastically (GenerativeDistribution).
 * A generator also takes named Parameters which are themselves Values, which may have been generated by a Generator.
 */
public interface Generator<T> extends GraphicalModelNode<T> {

    String getName();

    /**
     * @return a value generated by this generator.
     */
    Value<T> generate();

    String codeString();

    /**
     * @return the character symbol, for function '=' and for generative distribution '~'
     */
    char generatorCodeChar();

    @Override
    default List<GraphicalModelNode> getInputs() {
        return new ArrayList<>(getParams().values());
    }

    Map<String, Value> getParams();

    default String getInferenceStatement(Value value, Narrative narrative) {

        StringBuilder builder = new StringBuilder();

        builder.append("P(");

        String name = narrative.getId(value, false);

        builder.append(name);
        Map<String, Value> params = getParams();

        List<ParameterInfo> parameterInfos = getParameterInfo(0);
        int count = 0;
        for (ParameterInfo parameterInfo : parameterInfos) {
            Value v = params.get(parameterInfo.name());
            if (v != null && v.isRandom()) {

                if (count == 0) {
                    builder.append(" | ");
                } else {
                    builder.append(", ");
                }

                if (v.isAnonymous()) {
                    builder.append(parameterInfo.name());
                } else {

                    name = narrative.getId(v, false);
                    builder.append(name);
                }
                count += 1;
            }
        }
        builder.append(")");

        return builder.toString();
    }

    default String getInferenceNarrative(Value value, boolean unique, Narrative narrative) {

        String narrativeName = getNarrativeName();

        GeneratorInfo info = getGeneratorInfo(this.getClass());
        String citationString = narrative.cite(getCitation());

        String verbClause = info != null ? info.verbClause() : "comes from";
        StringBuilder builder = new StringBuilder();
        builder.append(NarrativeUtils.getValueClause(value, unique, narrative));
        builder.append(" ");
        builder.append(verbClause);
        builder.append(" ");
        if (!(this instanceof ExpressionNode)) {
            if (this instanceof DeterministicFunction) {
                builder.append(NarrativeUtils.getDefiniteArticle(narrativeName, true));
            } else {
                builder.append(NarrativeUtils.getIndefiniteArticle(narrativeName, true));
            }
        }
        builder.append(" ");
        builder.append(narrativeName);
        if (citationString != null && citationString != "") {
            builder.append(" ");
            builder.append(citationString);
        }

        Map<String, Value> params = getParams();
        String currentVerb = "";
        List<ParameterInfo> parameterInfos = getParameterInfo(0);
        int count = 0;
        for (ParameterInfo parameterInfo : parameterInfos) {
            Value v = params.get(parameterInfo.name());
            if (v != null) {
                if (count == 0) builder.append(" ");
                if (count > 0) {
                    if (count == params.size() - 1) {
                        builder.append(" and ");
                    } else {
                        builder.append(", ");
                    }
                }
                if (!parameterInfo.verb().equals(currentVerb)) {
                    currentVerb = parameterInfo.verb();
                    builder.append(currentVerb);
                    builder.append(" ");
                }
                builder.append(NarrativeUtils.getValueClause(v, false, true, false, this, narrative));
                count += 1;
            }
        }
        builder.append(".");
        return builder.toString();
    }

    /**
     * Get the name of the type of object generated by this generator.
     * @return
     */
    default String getTypeName() {
        return getReturnType(this.getClass()).getSimpleName();
    }

    default void setParam(String paramName, Value<?> value) {

        String methodName = "set" + Character.toUpperCase(paramName.charAt(0)) + paramName.substring(1);

        try {
            Method method = getClass().getMethod(methodName, value.value().getClass());

            method.invoke(this, value.value());
        } catch (NoSuchMethodException e) {

            Method[] methods = getClass().getMethods();
            for (Method method : methods) {
                if (method.getName().equals(methodName)) {
                    try {
                        method.invoke(this, value.value());
                        break;
                    } catch (InvocationTargetException | IllegalAccessException ignored) {
                    }
                }
            }
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        } catch (InvocationTargetException e) {
            e.printStackTrace();
        }
    }

    default void setInput(String paramName, Value<?> value) {
        setParam(paramName, value);
        value.addOutput(this);
    }

    default void setInputs(Map<String, Value<?>> params) {
        params.forEach(this::setInput);
    }

    default String getParamName(Value value) {
        Map<String, Value> params = getParams();
        for (String key : params.keySet()) {
            if (params.get(key) == value) return key;
        }
        return null;
    }

    /**
     * @return true if any of the parameters are random variables,
     * or are themselves that result of a function with random parameters as arguments.
     */
    default boolean hasRandomParameters() {
        for (Map.Entry<String, Value> entry : getParams().entrySet()) {

            Value<?> v = entry.getValue();

            if (v == null) {
                throw new RuntimeException("Unexpected null value for param " + entry.getKey() + " in generator " + getName());
            }

            if (v.isRandom()) return true;
        }
        return false;
    }

    default String getParamName(int paramIndex, int constructorIndex) {
        return getParameterInfo(constructorIndex).get(paramIndex).name();
    }

    @Deprecated
    default String getParamName(int paramIndex) {
        return getParamName(paramIndex, 0);
    }

    default List<ParameterInfo> getParameterInfo(int constructorIndex) {
        return getParameterInfo(this.getClass(), constructorIndex);
    }

    default Citation getCitation() {
        return getCitation(getClass());
    }

    default Class<?> getParamType(String name) {
        return getParams().get(name).getType();
    }

    default GeneratorInfo getInfo() {

        Class<?> classElement = getClass();

        Method[] methods = classElement.getMethods();

        for (Method method : methods) {
            for (Annotation annotation : method.getAnnotations()) {
                if (annotation instanceof GeneratorInfo) {
                    return (GeneratorInfo) annotation;
                }
            }
        }
        return null;
    }

    default String getRichDescription(int index) {

        List<ParameterInfo> pInfo = getParameterInfo(index);

        Map<String, Value> paramValues = getParams();

        StringBuilder html = new StringBuilder("<html><h3>");
        html.append(getName());
        if (this instanceof GenerativeDistribution) {
            html.append(" distribution");
        }
        html.append("</h3>");
        GeneratorInfo info = getInfo();
        if (info != null) {
            html.append("<p>").append(getInfo().description()).append("</p>");
        }
        if (pInfo.size() > 0) {
            html.append("<p>parameters: <ul>");
            for (ParameterInfo pi : pInfo) {
                html.append("<li>").append(pi.name()).append(" (").append(paramValues.get(pi.name())).append("); <font color=\"#808080\">").append(pi.description()).append("</font></li>");
            }
            html.append("</ul>");
        }

        Citation citation = getCitation();
        if (citation != null) {
            html.append("<h3>Reference</h3>");
            html.append(citation.value());
            String url = NarrativeUtils.getURL(citation);
            if (url.length() > 0) {
                html.append("<br><a href=\"" + url + "\">" + url + "</a><br>");
            }
        }

        html.append("</p></html>");
        return html.toString();
    }

    static List<ParameterInfo> getParameterInfo(Class<?> c, int constructorIndex) {
        return getParameterInfo(c.getConstructors()[constructorIndex]);
    }

    static Citation getCitation(Class<?> c) {
        Annotation[] annotations = c.getAnnotations();
        for (Annotation annotation : annotations) {
            if (annotation instanceof Citation) {
                return (Citation) annotation;
            }
        }
        return null;
    }


    static List<ParameterInfo> getParameterInfo(Constructor constructor) {
        ArrayList<ParameterInfo> pInfo = new ArrayList<>();

        Annotation[][] annotations = constructor.getParameterAnnotations();
        for (int i = 0; i < annotations.length; i++) {
            Annotation[] annotations1 = annotations[i];
            for (Annotation annotation : annotations1) {
                if (annotation instanceof ParameterInfo) {
                    pInfo.add((ParameterInfo) annotation);
                }
            }
        }

        return pInfo;
    }


    static Class<?>[] getParameterTypes(Class<? extends Generator> c, int constructorIndex) {
        return getParameterTypes(c.getConstructors()[constructorIndex]);
    }

    /**
     * @param constructor
     * @return an array of the generic types of arguments of the given constructor.
     */
    static Class[] getParameterTypes(Constructor constructor) {
        Type[] generics = constructor.getGenericParameterTypes();
        Class[] types = new Class[generics.length];
        for (int i = 0; i < generics.length; i++) {
            types[i] = lphy.reflection.Utils.getClass(generics[i]);
        }
        return types;
    }

    static List<Argument> getArguments(Class<?> c, int constructorIndex) {
        return getArguments(c.getConstructors()[constructorIndex]);
    }

    static List<Argument> getArguments(Constructor constructor) {

        List<Argument> arguments = new ArrayList<>();

        Annotation[][] annotations = constructor.getParameterAnnotations();
        Class<?>[] parameterTypes = getParameterTypes(constructor);

        // top for loop
        for (int i = 0; i < annotations.length; i++) {
            Annotation[] annotations1 = annotations[i];
            for (Annotation annotation : annotations1) {
                if (annotation instanceof ParameterInfo) {
                    arguments.add(new Argument(i, (ParameterInfo) annotation, parameterTypes[i]));
                }
            }
        }
        return arguments;
    }

    static String getGeneratorMarkdown(Class<? extends Generator> generatorClass) {

        GeneratorInfo generatorInfo = getGeneratorInfo(generatorClass);

        List<ParameterInfo> pInfo = getParameterInfo(generatorClass, 0);
        Class[] types = getParameterTypes(generatorClass, 0);

        StringBuilder md = new StringBuilder();

        StringBuilder signature = new StringBuilder();

        signature.append(Generator.getGeneratorName(generatorClass)).append("(");

        int count = 0;
        for (int i = 0; i < pInfo.size(); i++) {
            ParameterInfo pi = pInfo.get(i);
            if (count > 0) signature.append(", ");
            signature.append(new Text(types[i].getSimpleName())).append(" ").append(new BoldText(pi.name()));
            count += 1;
        }
        signature.append(")");

        md.append(new Heading(signature.toString(), 2)).append("\n\n");

        if (generatorInfo != null) md.append(generatorInfo.description()).append("\n\n");

        if (pInfo.size() > 0) {
            md.append(new Heading("Parameters", 3)).append("\n\n");
            List<Object> paramText = new ArrayList<>();

            for (int i = 0; i < pInfo.size(); i++) {
                ParameterInfo pi = pInfo.get(i);
                paramText.add(new Text(types[i].getSimpleName() + " " + new BoldText(pi.name()) + " - " + pi.description()));
            }
            md.append(new UnorderedList<>(paramText));
        }
        md.append("\n\n");

        md.append(new Heading("Return type", 3)).append("\n\n");

        List<String> returnType = Collections.singletonList(getReturnType(generatorClass).getSimpleName());
        md.append(new UnorderedList<>(returnType)).append("\n\n");

        Citation citation = getCitation(generatorClass);
        if (citation != null) {
            md.append(new Heading("Reference", 3)).append("\n\n");
            md.append(citation.value());
            if (citation.DOI().length() > 0) {
                String url = citation.DOI();
                if (!url.startsWith("http")) {
                    url = "http://doi.org/" + url;
                }
                md.append(new Link(url, url));
            }
        }
        return md.toString();
    }

    static String getGeneratorHtml(Class<? extends Generator> generatorClass) {
        GeneratorInfo generatorInfo = getGeneratorInfo(generatorClass);

        List<ParameterInfo> pInfo = getParameterInfo(generatorClass, 0);
        Class[] types = getParameterTypes(generatorClass, 0);

        // parameters
        StringBuilder signature = new StringBuilder();
        signature.append(Generator.getGeneratorName(generatorClass)).append("(");

        int count = 0;
        for (int i = 0; i < pInfo.size(); i++) {
            ParameterInfo pi = pInfo.get(i);
            if (count > 0) signature.append(", ");
            signature //.append(types[i].getSimpleName()).append(" ")
                    .append("<i>").append(pi.name()).append("</i>");
            count += 1;
        }
        signature.append(")");

        // main content
        StringBuilder html = new StringBuilder("<html><h2>");
        html.append(signature);
        html.append("</h2>");

        if (generatorInfo != null) html.append("<p>").append(generatorInfo.description()).append("</p>");

        if (pInfo.size() > 0) {
            html.append("<h3>Parameters:</h3>").append("<ul>");
//            int count = 0;
            for (int i = 0; i < pInfo.size(); i++) {
                ParameterInfo pi = pInfo.get(i);
                html.append("<li>").append(types[i].getSimpleName()).
                        append(" <b>").append(pi.name()).append("</b>")
                        .append(" - <font color=\"#808080\">")
                        .append(pi.description()).append("</font></li>");

//                if (count > 0) signature.append(", ");
//                signature.append(new Text(types[i].getSimpleName())).append(" ").append(new BoldText(pi.name()));
//                count += 1;
            }
//            signature.append(")");
//            html.append(new Heading(signature.toString(), 2)).append("\n\n");
            html.append("</ul>");
        }

        List<String> returnType = Collections.singletonList(getReturnType(generatorClass).getSimpleName());
        if (returnType.size() > 0) {
            html.append("<h3>Return type:</h3>").append("<ul>");
            for (String itm : returnType)
                html.append("<li>").append(itm).append("</li>");
            html.append("</ul>");
        }

        Citation citation = getCitation(generatorClass);
        if (citation != null) {
            html.append("<h3>Reference</h3>");
            html.append(citation.value());
            String url = NarrativeUtils.getURL(citation);
            if (url.length() > 0)
                html.append("<br><a href=\"").append(url).append("\">").append(url).append("</a><br>");
        }

        html.append("</p></html>");
        return html.toString();
    }

    static List<ParameterInfo> getAllParameterInfo(Class c) {
        ArrayList<ParameterInfo> pInfo = new ArrayList<>();
        for (Constructor constructor : c.getConstructors()) {
            pInfo.addAll(getParameterInfo(constructor));
        }
        return pInfo;
    }

    static String getSignature(Class<?> aClass) {

        List<ParameterInfo> pInfo = Generator.getParameterInfo(aClass, 0);

        StringBuilder builder = new StringBuilder();
        builder.append(getGeneratorName(aClass));
        builder.append("(");
        if (pInfo.size() > 0) {
            builder.append(pInfo.get(0).name());
            for (int i = 1; i < pInfo.size(); i++) {
                builder.append(", ");
                builder.append(pInfo.get(i).name());
            }
        }
        builder.append(")");
        return builder.toString();
    }

    static String getGeneratorName(Class<?> c) {
        GeneratorInfo ginfo = getGeneratorInfo(c);
        if (ginfo != null) return ginfo.name();
        return c.getSimpleName();
    }

    static String getGeneratorDescription(Class<?> c) {
        GeneratorInfo ginfo = getGeneratorInfo(c);
        if (ginfo != null) return ginfo.description();
        return "";
    }

    static GeneratorInfo getGeneratorInfo(Class<?> c) {

        Method[] methods = c.getMethods();
        for (Method method : methods) {
            for (Annotation annotation : method.getAnnotations()) {
                if (annotation instanceof GeneratorInfo) {
                    return (GeneratorInfo) annotation;
                }
            }
        }
        return null;
    }

    static String getArgumentCodeString(Map.Entry<String, Value> entry) {
        return getArgumentCodeString(entry.getKey(), entry.getValue());
    }

    static String getArgumentCodeString(String name, Value value) {
        String prefix = "";
        if (!Utils.isInteger(name)) {
            prefix = name + "=";
        }

        if (value == null) {
            throw new RuntimeException("Value of " + name + " is null!");
        }

        if (value.isAnonymous()) return prefix + value.codeString();
        return prefix + value.getId();
    }

    /**
     * @param arguments the arguments of the Generator
     * @param initArgs the parallel array of initial arguments that match the arguments of the generator - may contain nulls where no name match
     * @param params the map of all params, may be more than non-null in initial arguments
     * @param lightweight
     * @return
     */
    static boolean matchingParameterTypes(List<Argument> arguments, Object[] initArgs, Map<String, Value> params, boolean lightweight) {

        int count = 0;
        for (int i = 0; i < arguments.size(); i++) {
            Argument argumentInfo = arguments.get(i);
            Object arg = initArgs[i];

            if (arg != null) {
                Class parameterType = argumentInfo.type;
                Class valueType = lightweight ? arg.getClass() : ((Value) arg).value().getClass();

                if (!parameterType.isAssignableFrom(valueType)) return false;
                count += 1;
            } else {
                if (!argumentInfo.optional) return false;
            }
        }
        return params == null || count == params.size();
    }

    static Map<String, Value> convertArgumentsToParameterMap(List<Argument> argumentInfos, Object[] initArgs) {
        Map<String, Value> params = new TreeMap<>();
        for (int i = 0; i < argumentInfos.size(); i++) {
            Argument argumentInfo = argumentInfos.get(i);
            Value value = (Value) initArgs[i];

            if (value != null) params.put(argumentInfo.name, value);
        }
        return params;
    }

    static Class<?> getReturnType(Class<?> genClass) {
        Method[] methods = genClass.getMethods();

        for (Method method : methods) {
            GeneratorInfo generatorInfo = method.getAnnotation(GeneratorInfo.class);
            if (generatorInfo != null) {
                return lphy.reflection.Utils.getGenericReturnType(method);
            }
        }
        if (GenerativeDistribution.class.isAssignableFrom(genClass)) {
            try {
                Method method = genClass.getMethod("sample");
                return lphy.reflection.Utils.getGenericReturnType(method);
            } catch (NoSuchMethodException e) {
                e.printStackTrace();
            }
        } else if (DeterministicFunction.class.isAssignableFrom(genClass)) {
            {
                try {
                    Method method = genClass.getMethod("apply");
                    return lphy.reflection.Utils.getGenericReturnType(method);
                } catch (NoSuchMethodException e) {
                    e.printStackTrace();
                }
            }
        }
        return Object.class;
    }

    /**
     * @param value
     * @return the narrative name for the given value, being a parameter of this generator.
     */
    default String getNarrativeName(Value value) {
        return getNarrativeName(getParamName(value));
    }

    /**
     * @param paramName the parameter name
     * @return the narrative name for the given parameter name.
     */
    default String getNarrativeName(String paramName) {
        List<ParameterInfo> parameterInfos = getParameterInfo(0);
        for (ParameterInfo parameterInfo : parameterInfos) {
            if (parameterInfo.name().equals(paramName)) {
                if (parameterInfo.suppressNameInNarrative()) return "";
                if (parameterInfo.narrativeName().length() > 0) {
                    return parameterInfo.narrativeName();
                }
            }
        }
        return paramName;
    }

    /**
     * @return the narrative name of this generator.
     */
    default String getNarrativeName() {
        GeneratorInfo generatorInfo = getGeneratorInfo(this.getClass());
        if (generatorInfo != null) {
            if (generatorInfo.narrativeName().length() > 0) return generatorInfo.narrativeName();
        }
        return getName();
    }
}
