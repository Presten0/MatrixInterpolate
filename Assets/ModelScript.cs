using System;
using UnityEditor;
using UnityEngine;
using Vectors;

[ExecuteAlways]
[RequireComponent(typeof(VectorRenderer))]
public class ModelScript : MonoBehaviour {
    
    [NonSerialized] 
    private VectorRenderer vectors;

    [SerializeField, Range(0f, 1f)] 
    private float Time;
    
    [SerializeField,HideInInspector]internal Matrix4x4 matrixA = Matrix4x4.identity;
    [SerializeField,HideInInspector]internal Matrix4x4 matrixB = Matrix4x4.identity;
    [SerializeField,HideInInspector]internal Matrix4x4 matrixI = Matrix4x4.identity;

    [SerializeField] 
    private bool interpolateTranslation = true;
    [SerializeField] 
    private bool interpolateRotation = true;
    [SerializeField] 
    private bool interpolateScale = true;

    void OnEnable() {
        vectors = GetComponent<VectorRenderer>();
    }

    void Update()
    {
        //matrixA = Matrix4x4.identity;
        //matrixB = Matrix4x4.identity;

        if (interpolateTranslation) //Check if user wants to interpolate translation
        {
            var interpolatedTranslation = (1.0f - Time) * GetTranslationVector(matrixA) + Time * GetTranslationVector(matrixB);
            SetTranslationVector(ref matrixI, interpolatedTranslation);
        }
        else
        {
            SetTranslationVector(ref matrixI, GetTranslationVector(matrixA));
        }

        if (interpolateRotation) //Check if user wants to interpolate rotation
        {
            var rotationA = GetRotation(matrixA);
            var rotationB = GetRotation(matrixB);
            var rotationI = GetInterpolatedQuaternion(rotationA, rotationB);
            var rotationMatrixI = ConvertQuaternion(rotationI);
            matrixI = SetRotation( matrixI, rotationMatrixI);
        }
        else
        {
            matrixI = SetRotation(matrixI, ConvertQuaternion(GetRotation(matrixA)));
        }

        if (interpolateScale) //Check if user wants to interpolate scale
        {
            SetScale(ref matrixI, GetInterpolatedScaleMagnitudes(matrixA, matrixB));
        }
        else
        {
            SetScale(ref matrixI, new Vector3(1,1,1));
        }
        
        // Variables for vector addition used for the visualization of the unit vectors per matrix
        
        Vector3 AX = matrixA.GetColumn(0); //GetColumn returns a vector representing the column with index, example index:0 returns (m00, m10, m20)
        Vector3 AY = matrixA.GetColumn(1); 
        Vector3 AZ = matrixA.GetColumn(2);
        
        Vector3 BX = matrixB.GetColumn(0);
        Vector3 BY = matrixB.GetColumn(1);
        Vector3 BZ = matrixB.GetColumn(2);
        
        Vector3 IX = matrixI.GetColumn(0);
        Vector3 IY = matrixI.GetColumn(1);
        Vector3 IZ = matrixI.GetColumn(2);
        
        using (vectors.Begin()) {
            
            //Visualize matrix A
            vectors.Draw(GetTranslationVector(matrixA), GetTranslationVector(matrixA)+AX, Color.red);
            vectors.Draw(GetTranslationVector(matrixA), GetTranslationVector(matrixA)+AY, Color.green);
            vectors.Draw(GetTranslationVector(matrixA), GetTranslationVector(matrixA)+AZ, Color.blue);
            
            //Visualize matrix B
            vectors.Draw(GetTranslationVector(matrixB), GetTranslationVector(matrixB)+BX, Color.red);
            vectors.Draw(GetTranslationVector(matrixB), GetTranslationVector(matrixB)+BY, Color.green);
            vectors.Draw(GetTranslationVector(matrixB), GetTranslationVector(matrixB)+BZ, Color.blue);
            
            //Visualize box from interpolated matrix
            vectors.Draw(GetTranslationVector(matrixI), GetTranslationVector(matrixI)+IX, Color.red);
            vectors.Draw(GetTranslationVector(matrixI)+IY, GetTranslationVector(matrixI)+IX+IY, Color.red);
            vectors.Draw(GetTranslationVector(matrixI)+IZ, GetTranslationVector(matrixI)+IX+IZ, Color.red);
            vectors.Draw(GetTranslationVector(matrixI)+IY+IZ, GetTranslationVector(matrixI)+IX+IY+IZ, Color.red);
            vectors.Draw(GetTranslationVector(matrixI), GetTranslationVector(matrixI)+IY, Color.green);
            vectors.Draw(GetTranslationVector(matrixI)+IX, GetTranslationVector(matrixI)+IY+IX, Color.green);
            vectors.Draw(GetTranslationVector(matrixI)+IZ, GetTranslationVector(matrixI)+IY+IZ, Color.green);
            vectors.Draw(GetTranslationVector(matrixI)+IX+IZ, GetTranslationVector(matrixI)+IY+IX+IZ, Color.green);
            vectors.Draw(GetTranslationVector(matrixI), GetTranslationVector(matrixI)+IZ, Color.blue);
            vectors.Draw(GetTranslationVector(matrixI)+IX, GetTranslationVector(matrixI)+IZ+IX, Color.blue);
            vectors.Draw(GetTranslationVector(matrixI)+IY, GetTranslationVector(matrixI)+IZ+IY, Color.blue);
            vectors.Draw(GetTranslationVector(matrixI)+IX+IY, GetTranslationVector(matrixI)+IZ+IX+IY, Color.blue);
        }
    }

    public Vector3 GetTranslationVector(Matrix4x4 matrix)
    {
        return matrix.GetColumn(3); //Returns Vector3(m03, m13, m23)
    }
    
    public void SetTranslationVector(ref Matrix4x4 matrix, Vector3 translationVector) //Sets 
    {
        matrix.m03 = translationVector.x;
        matrix.m13 = translationVector.y;
        matrix.m23 = translationVector.z;
    }

    public Vector3 GetInterpolatedScaleMagnitudes(Matrix4x4 matrixA, Matrix4x4 matrixB)
    {
        var scaleAX = MathF.Sqrt(matrixA.m00 * matrixA.m00 + matrixA.m10 * matrixA.m10 + matrixA.m20 * matrixA.m20);
        var scaleAY = MathF.Sqrt(matrixA.m01 * matrixA.m01 + matrixA.m11 * matrixA.m11 + matrixA.m21 * matrixA.m21);
        var scaleAZ = MathF.Sqrt(matrixA.m02 * matrixA.m02 + matrixA.m12 * matrixA.m12 + matrixA.m22 * matrixA.m22);
        
        var scaleBX = MathF.Sqrt(matrixB.m00 * matrixB.m00 + matrixB.m01 * matrixB.m01 + matrixB.m02 * matrixB.m02);
        var scaleBY = MathF.Sqrt(matrixB.m01 * matrixB.m01 + matrixB.m11 * matrixB.m11 + matrixB.m21 * matrixB.m21);
        var scaleBZ = MathF.Sqrt(matrixB.m02 * matrixB.m02 + matrixB.m12 * matrixB.m12 + matrixB.m22 * matrixB.m22);

        var scaleIX = (1.0f - Time) * scaleAX + Time * scaleBX;
        var scaleIY = (1.0f - Time) * scaleAY + Time * scaleBY;
        var scaleIZ = (1.0f - Time) * scaleAZ + Time * scaleBZ;
        
        return new Vector3(scaleIX, scaleIY, scaleIZ); //Scaling vector
    }

    public Vector3 GetScale(Matrix4x4 matrix)
    {
        var scaleX = MathF.Sqrt(matrix.m00 * matrix.m00 + matrix.m10 * matrix.m10 + matrix.m20 * matrix.m20);
        var scaleY = MathF.Sqrt(matrix.m01 * matrix.m01 + matrix.m11 * matrix.m11 + matrix.m21 * matrix.m21);
        var scaleZ = MathF.Sqrt(matrix.m02 * matrix.m02 + matrix.m12 * matrix.m12 + matrix.m22 * matrix.m22);

        return new Vector3(scaleX, scaleY, scaleZ); //Scaling vector
    }

    public void SetScale(ref Matrix4x4 matrix, Vector3 scaleVector)
    {
        matrix.SetColumn(0,matrix.GetColumn(0)*scaleVector.x);
        matrix.SetColumn(1,matrix.GetColumn(1)*scaleVector.y);
        matrix.SetColumn(2,matrix.GetColumn(2)*scaleVector.z);
    }

    public Quaternion GetRotation(Matrix4x4 matrix)
    {
        var matrixColumnX = matrix.GetColumn(0);
        var matrixColumnY = matrix.GetColumn(1);
        var matrixColumnZ = matrix.GetColumn(2);

        matrix.SetColumn(0, matrixColumnX.normalized);
        matrix.SetColumn(1, matrixColumnY.normalized);
        matrix.SetColumn(2, matrixColumnZ.normalized);

        float fourWSquaredMinus1 = matrix.m00 + matrix.m11 + matrix.m22;
        float fourXSquaredMinus1 = matrix.m00 - matrix.m11 - matrix.m22;
        float fourYSquaredMinus1 = matrix.m11 - matrix.m00 - matrix.m22;
        float fourZSquaredMinus1 = matrix.m22 - matrix.m00 - matrix.m11;

        float w = 0, x = 0, y = 0, z = 0;
        
        int biggestIndex = 0;
        float fourBiggestSquaredMinus1 = fourWSquaredMinus1;
        
        if (fourXSquaredMinus1 > fourBiggestSquaredMinus1)
        {
            fourBiggestSquaredMinus1 = fourXSquaredMinus1;
            biggestIndex = 1;
        }

        if (fourYSquaredMinus1 > fourBiggestSquaredMinus1)
        {
            fourBiggestSquaredMinus1 = fourYSquaredMinus1;
            biggestIndex = 2;
        }

        if (fourZSquaredMinus1 > fourBiggestSquaredMinus1)
        {
            fourBiggestSquaredMinus1 = fourZSquaredMinus1;
            biggestIndex = 3;
        }

        // Perform square root and division
        float biggestVal = MathF.Sqrt(fourBiggestSquaredMinus1 + 1.0f) * 0.5f;
        float mult = 0.25f / biggestVal;
        // Apply table to compute quaternion values
        switch (biggestIndex)
        {
            case 0:
                w = biggestVal;
                x = (matrix.m21 - matrix.m12) * mult;
                y = (matrix.m02 - matrix.m20) * mult;
                z = (matrix.m10 - matrix.m01) * mult;
                break;

            case 1:
                x = biggestVal;
                w = (matrix.m21 - matrix.m12) * mult;
                y = (matrix.m10 - matrix.m01) * mult;
                z = (matrix.m02 - matrix.m20) * mult;
                break;

            case 2:
                y = biggestVal;
                w = (matrix.m02 - matrix.m20) * mult;
                x = (matrix.m10 - matrix.m01) * mult;
                z = (matrix.m21 - matrix.m12) * mult;
                break;

            case 3:
                z = biggestVal;
                w = (matrix.m10 - matrix.m01) * mult;
                x = (matrix.m02 - matrix.m20) * mult;
                y = (matrix.m21 - matrix.m12) * mult;
                break;

        }

        return new Quaternion(x, y, z, w);

    }


    public float ProjectPairAngle(Vector3 vectorPrim, Vector3 vectorUnit, Vector3 normal)
    {
        var normalMag = MathF.Sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z); //Get magnitude of normal of plane
        var projPrimOnNormal = ((Vector3.Dot(vectorPrim, normal) / (normalMag * normalMag))) * normal; //Project the vector onto the normal
        var projUnitOnNormal = ((Vector3.Dot(vectorUnit, normal) / (normalMag * normalMag))) * normal; //Project the vector onto the normal

        var projPrim = vectorPrim - projPrimOnNormal; //Vector-subtraction to get primVector parallel to plane (projection onto plane)
        var projUnit = vectorUnit - projUnitOnNormal; //Vector-subtraction to get unitVector parallel to plane (projection onto plane)
        
        var scalarProd = Vector3.Dot(projUnit, projPrim); //Dot-product between the projected vectors
        
        //Get magnitude of each vector projected onto plane
        var projprimMag = MathF.Sqrt((projPrim.x * projPrim.x) + (projPrim.y * projPrim.y) +
                                  (projPrim.z * projPrim.z));
        var projUnitMag = MathF.Sqrt((projUnit.x * projUnit.x) + (projUnit.y * projUnit.y) +
                                  (projUnit.z * projUnit.z));
        
        //Formula to convert dot-product to angle
        return MathF.Acos(scalarProd / (projprimMag * projUnitMag));
        
    }

    public Quaternion GetInterpolatedQuaternion(Quaternion quatA, Quaternion quatB)
    {
        // Variables for components
        float ax = quatA.x;
        float ay = quatA.y;
        float az = quatA.z;
        float aw = quatA.w;
        
        float bx = quatB.x;
        float by = quatB.y;
        float bz = quatB.z;
        float bw = quatB.w;

        float t = Time; // Interpolation variable

        float cx,cy,cz,cw; // Output components
        
        float cosOmega = aw*bw + ax*bx + ay*by + az*bz; // Get the cosine
        
        if (cosOmega < 0.0f) { // If the cosine is negative, negate components of one of the quaternions
            bw = -bw;
            bx = -bx;
            by = -by;
            bz = -bz;
            cosOmega = -cosOmega;
        }

        float k0, k1;
        if (cosOmega > 0.9999f) { // To not be able to divide by zero

            //linear interpolation
            k0 = 1.0f-t;
            k1 = t;

        } else {

            float sinOmega = Mathf.Sqrt(1.0f - cosOmega*cosOmega); // Compute the sin of the angle using the trig identity sin^2(omega) + cos^2(omega) = 1
            
            float omega = Mathf.Atan2(sinOmega, cosOmega); // Compute the angle from its sine and cosine

            float oneOverSinOmega = 1.0f / sinOmega; // Compute inverse of denominator

            // Compute interpolation parameters
            k0 = Mathf.Sin((1.0f - t) * omega) * oneOverSinOmega;
            k1 = Mathf.Sin(t * omega) * oneOverSinOmega;
        }
        
        // Interpolation
        cx = ax*k0 + bx*k1;
        cy = ay*k0 + by*k1;
        cz = az*k0 + bz*k1;
        cw = aw*k0 + bw*k1;

        return new Quaternion(cx, cy, cz, cw);
    }

    public Matrix4x4 ConvertQuaternion(Quaternion quat)
    {
        var newMatrix = Matrix4x4.identity;
        
        //Column 1
        newMatrix.m00 = 1 - 2 * (quat.y * quat.y) - 2 * (quat.z * quat.z);
        newMatrix.m10 = 2 * (quat.x * quat.y) - 2 * (quat.w * quat.z);
        newMatrix.m20 = 2 * (quat.x * quat.z) + 2 * (quat.w * quat.y);
        //Column 2
        newMatrix.m01 = 2 * (quat.x * quat.y) + 2 * (quat.w * quat.z);
        newMatrix.m11 = 1 - 2 * (quat.x * quat.x) - 2 * (quat.z * quat.z);
        newMatrix.m21 = 2 * (quat.y * quat.z) - 2 * (quat.w * quat.x);
        //Column 3
        newMatrix.m02 = 2 * (quat.x * quat.z) - 2 * (quat.w * quat.y);
        newMatrix.m12 = 2 * (quat.y * quat.z) + 2 * (quat.w * quat.x);
        newMatrix.m22 = 1 - 2 * (quat.x * quat.x) - 2 * (quat.y * quat.y);

        return newMatrix.transpose;
    }

    public Matrix4x4 SetRotation(Matrix4x4 matrix, Matrix4x4 rotationMatrix)
    {
        matrix.m00 = rotationMatrix.m00;
        matrix.m01 = rotationMatrix.m01;
        matrix.m02 = rotationMatrix.m02;
        matrix.m10 = rotationMatrix.m10;
        matrix.m11 = rotationMatrix.m11;
        matrix.m12 = rotationMatrix.m12;
        matrix.m20 = rotationMatrix.m20;
        matrix.m21 = rotationMatrix.m21;
        matrix.m22 = rotationMatrix.m22;

        return matrix;
    }

    public Vector3 GetNormalizedVector(Vector3 vector)
    {
        var magnitude = MathF.Sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
        return vector / magnitude;
    }
}

[CustomEditor(typeof(ModelScript))]
public class ModelGUI : Editor {
    void OnSceneGUI() {
        var ex = target as ModelScript;
        if (ex == null) return;

        EditorGUI.BeginChangeCheck();
        var targetA = Handles.PositionHandle(ex.GetTranslationVector(ex.matrixA), ex.GetRotation(ex.matrixA)); //Handle at translation-vector
        var targetB = Handles.PositionHandle(ex.GetTranslationVector(ex.matrixB), ex.GetRotation(ex.matrixB)); //Handle at translation-vector
        
        Debug.Log(ex.GetRotation(ex.matrixB));
        
        if (EditorGUI.EndChangeCheck()) {
            Undo.RecordObject(target, "Vector Positions");

            var copyA = ex.matrixA;
            var copyB = ex.matrixB;
            
            ex.SetTranslationVector(ref copyA, targetA); //Set translation-vector to handles position
            ex.SetTranslationVector(ref copyB, targetB); //Set translation-vector to handles position

            ex.matrixA = copyA;
            ex.matrixB = copyB;

            EditorUtility.SetDirty(target);
        }
        
        EditorGUI.BeginChangeCheck();
        var c = Handles.RotationHandle(ex.GetRotation(ex.matrixA), ex.GetTranslationVector(ex.matrixA));
        var d = Handles.RotationHandle(ex.GetRotation(ex.matrixB), ex.GetTranslationVector(ex.matrixB));

        if (EditorGUI.EndChangeCheck()) {
            Undo.RecordObject(target, "Vector Positions");
            ex.matrixA = ex.SetRotation(ex.matrixA, ex.ConvertQuaternion(c));
            ex.matrixB = ex.SetRotation(ex.matrixB, ex.ConvertQuaternion(d));

            EditorUtility.SetDirty(target);
        }
    }
    public override void OnInspectorGUI()
    {
        base.OnInspectorGUI();
        
        var ex = target as ModelScript;
        if (ex == null) return;

        EditorGUI.BeginChangeCheck();

        EditorGUILayout.BeginHorizontal();
        EditorGUILayout.PrefixLabel("Matrix A");
        EditorGUILayout.BeginVertical();

        var resultA = Matrix4x4.identity;
        for (var i = 0; i < 4; i++)
        {
            EditorGUILayout.BeginHorizontal();
            for (var j = 0; j < 4; j++)
            {
                resultA[i, j] = EditorGUILayout.FloatField((ex.matrixA[i, j]));
            }
            EditorGUILayout.EndHorizontal();
        }
        
        EditorGUILayout.EndVertical();
        EditorGUILayout.EndHorizontal();
        
        EditorGUILayout.Space();
        
        EditorGUILayout.BeginHorizontal();
        EditorGUILayout.PrefixLabel("Matrix B");
        EditorGUILayout.BeginVertical();

        var resultB = Matrix4x4.identity;
        for (var k = 0; k < 4; k++)
        {
            EditorGUILayout.BeginHorizontal();
            for (var l = 0; l < 4; l++)
            {
                resultB[k, l] = EditorGUILayout.FloatField((ex.matrixB[k, l]));
            }
            EditorGUILayout.EndHorizontal();
        }
        
        EditorGUILayout.EndVertical();
        EditorGUILayout.EndHorizontal();
        
        EditorGUILayout.Space();
        
        EditorGUILayout.BeginHorizontal();
        EditorGUILayout.PrefixLabel("Matrix I");
        EditorGUILayout.BeginVertical();

        var resultI = Matrix4x4.identity;
        for (var o = 0; o < 4; o++)
        {
            EditorGUILayout.BeginHorizontal();
            for (var p = 0; p < 4; p++)
            {
                resultI[o, p] = EditorGUILayout.FloatField((ex.matrixI[o, p]));
            }
            EditorGUILayout.EndHorizontal();
        }
        
        EditorGUILayout.EndVertical();
        EditorGUILayout.EndHorizontal();
        
        if (EditorGUI.EndChangeCheck())
        {
            Undo.RecordObject(ex, "Changed matrix");
            ex.matrixA = resultA;
            ex.matrixB = resultB;
            EditorUtility.SetDirty(ex);
        }
    }
}


